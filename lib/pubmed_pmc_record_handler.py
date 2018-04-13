
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for handling data from biopython files for pubmed abstracts and pmc full-text articles.
Author: Tilia Ellendorff

From pmc full-text articles in-text tags are removed. These are checked towards a list of known relevant tags.

nxml is stored for pmc papers containing full text as well as those not containing fulltext if nxml_dump option is used.


"""

import codecs
import copy
import sys
import pickle
import os
from Bio import Entrez
import requests
import re
import argparse
import logging
import traceback
import unicodedata
from lxml import etree as et
import pandas as pd

from abbrev_matcher import AbbrevMatcher


class PMID_PMC_Mapper(object):
    '''
    ID mapping object: loads/saves/updates and applies dictionary mapping PMIDs to PMCs and vice-versa.
    ::data_path: path to where the dictionary is stored
    '''
    def __init__(self, data_path=None):
        self.data_path = data_path
        self.map_count = 0  # count overall number of mappings performed
        self.download_count = 0  # count number of downloads performed
        self.new_data = dict()  # dictionary of newly downloaded records

        if self.data_path is None:
            self.store_data = False
            self.data_format = None
            self.data = dict()
            print('NO LOCATION GIVEN FOR PMID PMCID MAPPING DATA. DATA WILL NOT BE STORED.')
        else:
            self.store_data = True
            if self.data_path.endswith('.pkl'):
                self.data_format = 'pickle'
                self.header = None
            elif self.data_path.endswith('.tsv'):
                self.data_format = 'table'
                self.header = ['in_id', 'out_id']
            else:
                print('Error: unrecognized format of id mapper data, please provide .pkl or .tsv path')
                raise Exception

        if self.data_format is not None:
            try:
                self.data = self._load_data(self.data_path)
            except IOError:
                print('NO STORED PMID PMCID MAPPING DATA FOUND UNDER THE GIVEN LOCATION. CREATING NEW ONE.')
                self.data = {}
            
    def _load_data(self, data_path):
        if self.data_format == 'pickle':
            data = pickle.load(open(data_path, 'rb'))
            return data
        if self.data_format == 'table':
            dict_df = pd.read_csv(self.data_path, sep='\t', dtype=str, names=self.header, index_col=None)
            dict_df.fillna('None', inplace=True)
            data = dict(zip(dict_df['in_id'], dict_df['out_id']))
            print('ID mapping dict loaded from table.')
            return data

    def _table_add_new(self):
        new_data_items = self.new_data.items()
        dict_df = pd.DataFrame(list(new_data_items), columns=self.header)
        dict_df.reindex_axis(self.header, axis=1)
        dict_df.to_csv(self.data_path, mode='a', sep='\t', columns=self.header, encoding='utf-8', header=False, index=False)
        self.new_data = dict()

    def save_data(self):
        if not self.store_data:
            print('Error: Data cannot be stored, no path is given.')
            pass
        elif self.data_format == 'pickle':
            pickle.dump(self.data_path, open(self.pickle_loc, 'wb'))
        else:
            mapper_items = self.data.items()
            dict_df = pd.DataFrame(list(mapper_items), columns=['in_id', 'out_id'])
            dict_df.sort_values(['in_id'], inplace=True)
            print('Shape of Output Dataframe', dict_df.shape)
            dict_df.to_csv(self.data_path, sep='\t', encoding='utf-8', index=False)
        print('Id mapping data has been saved under', self.data_path)
       
    def map_ids(self, available_id, entrez_email=None):
        '''Maps PMC IDs to Pubmed ids (PMIDs) and vice versa.'''
        self.map_count += 1
        try:
            result = self.data[available_id]

            if result is 'None' or result is None:
                return None

            else:
                return str(result)
                
        except KeyError:
            rq = 'http://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=' + \
                str(available_id) + '&tool=Python&email=' + str(entrez_email)
            try:
                r = requests.get(rq)
                self.download_count += 1
            except Exception as e:
                print('Download Problem')
                self.save_data()  # save data if download error occurs
                logging.error(traceback.format_exc())
                raise Exception("Download Error")
            try:
                tree = et.fromstring(r.text)
            except et.XMLSyntaxError:
                print('XML parsing error in id mapper:', rq)
                return None
            except ValueError:
                print('Problem with id mapper xml processing:', r.text)
                return None
            
            #print(et.tostring(tree, pretty_print=True))
            
            for record in tree.iterfind('.//record'):
                pmid = record.get("pmid")
                pmcid = record.get("pmcid")
                
                self.data[pmid] = pmcid
                self.data[pmcid] = pmid

                # adding to new data to be written out regularly
                self.new_data[pmid] = pmcid
                self.new_data[pmcid] = pmid
                if self.download_count % 50 == 0:  # save data after every 50 downloads
                    self._table_add_new()  # adding rows with new data to the table

            try:
                return self.data[available_id]
            except KeyError:
                # in case the online conversion tool throws an error message
                print('Invalid Article ID Error.')
                return None


class IDCollection(object):
    '''
    Class for reading a collection of file ids stored in a text file (one line per id) and reading them from the dump dir.
    '''
    def __init__(self, dump_dir, type='pmid', id_file=False, id_mapper=None, entrez_email=None, debug=False):
        '''type: pmid or pmcid; pmid is default'''

        self.dumpfiles_list = None

        self.type = type
        if self.type not in ['pmid', 'pmcid']:
            print('Unknown collection type:', self.type)
            raise AttributeError
        self.dump_dir = dump_dir
        if id_file:
            self.ids = [id.rstrip('\n') for id in open(id_file, 'r').readlines()]
        else:
            print('Attention: NO IDS GIVEN FOR COLLECTION')
            self.ids = None

        self.num_ids = len(self.ids)

        self.entrez_email = entrez_email

        self.dumpfiles_list = self.dumpfiles(id_mapper=id_mapper, debug=debug)


    def dumpfiles(self, id_mapper=None, debug=False):
        '''Returns a list of dumpfile objects for the defined collection type.'''
        if self.dumpfiles_list:
            return self.dumpfiles_list
        else:
            dumpfiles = list()

            if self.ids:
                for doc_id in self.ids:
                    if self.type == 'pmid':
                        id_dump_file = PMID_dump_file(
                            doc_id,
                            self.dump_dir,
                            id_mapper=id_mapper,
                            entrez_email=self.entrez_email,
                            debug=debug
                                                    )
                    elif self.type == 'pmcid':
                        id_dump_file = PMC_dump_file(
                            doc_id,
                            self.dump_dir,
                            id_mapper=id_mapper,
                            entrez_email=self.entrez_email,
                            debug=debug
                                                    )
                    dumpfiles.append(id_dump_file)

                return dumpfiles

    def get_section_collection(
            self,
            sections_list,
            out_path=None,
            separate=False,
            titles=False,
            id_mapper=None,
            report=False,
            debug=False
                            ):
        '''
        Only use specific specified sections of text, e.g. CONCLUSION and RESULTS.
        So far only for pmid collection type.
        Args:
            separate:
                False = write text for the whole collection into one file;
                True = write text for each document into a seperate file (in this case 'out path' is used as directory)
        '''

        text_list = list()  # stores text of the whole collection

        if separate and out_path is not None:  # write section text for each document to a separate file
            if not os.path.isdir(out_path):
                os.makedirs(out_path)
            write_file = False
            write_dir = True
        elif not separate and out_path is not None:  # write section text for all documents to one joint file
            write_file = True
            write_dir = False
            out_file = open(out_path, 'w', encoding='utf-8')
        else:  # do not write anything to a file
            print('NO OUTPUT WILL BE WRITTEN')
            write_file = False
            write_dir = False


        counter_dict = {sec:0 for sec in sections_list}

        for df in self.dumpfiles_list:
            abs_sec_dict = df.get_abstract_sections()
            one_doc_list = list()
            for sec in sections_list:
                try:
                    sec_text = abs_sec_dict[sec]
                    counter_dict[sec] += 1
                    if titles:
                        sec_text = sec + ': ' + sec_text
                    if out_path != None and write_file == True:
                        out_file.write(sec_text + '\n')
                    one_doc_list.append(sec_text)
                except KeyError:
                    pass

                if separate is False and one_doc_list:
                    text_list.extend(one_doc_list)
                elif separate is True and one_doc_list:
                    text_list.append(one_doc_list)

                if write_dir and one_doc_list:  # write to a separate file
                    filepath = os.path.join(out_path, df.doc_id)
                    out_file = open(filepath, 'w', encoding='utf-8')
                    for sec_text in one_doc_list:
                        out_file.write(sec_text + '\n')
                    out_file.close()

        if report is True:
            print(counter_dict)

        if write_file is True:
            out_file.close()

        return text_list

    def abstract_text_count(self):
        '''
        Counts how many ids in the collection have an abstract text included; 
        only valid for pmid collection type.
        '''
        if not self.type == 'pmid':
            print('Only possible with pmid collection type')
            raise ValueError


        count_dict = {'True':0, 'False':0}
        for df in self.dumpfiles_list:
            if df.has_abstract_text() == True:
                count_dict['True'] += 1
            else:
                count_dict['False'] += 1


        print('Total ids in the collection:', self.num_ids)
        print('Ids with an abstract text:', count_dict['True'])
        print('Ids without an abstract text:', count_dict['False'])


class DumpFile(object):
    '''Abstract DumpFile object. Includes PubMed_dump_files and PMC_dump_files.'''

    def __init__(self, doc_id, dump_dir, file_extension, id_mapper=None, entrez_email=None, debug=False):
        self.doc_id = doc_id
        self.dump_dir = dump_dir
        self.entrez_email = entrez_email
        self.file_extension = file_extension
        self.debug = debug

        if self.file_extension == '.nxml':
            self.type = 'pmc'
        elif self.file_extension == '.pxml':
            self.type = 'pubmed'
        else:
            print('Unknown file type:', self.file_extension)
            raise AttributeError

        if id_mapper != None:
            if doc_id.startswith('PMC'):
                self.pmc_id = self.doc_id
                self.pmid = id_mapper.map_ids(self.doc_id)
            else:
                self.pmid = doc_id
                self.pmc_id = id_mapper.map_ids(self.doc_id)
        else:
            if doc_id.startswith('PMC'):
                self.pmc_id = doc_id
                self.pmid = None
            else:
                self.pmid = doc_id
                self.pmc_id = None

        if self.type == 'pmc' and self.pmc_id == None:
            print('No PMC ID could be found for PMID', self.pmid)
            raise ValueError

        if self.type == 'pubmed':
            self.file_name = self.pmid + self.file_extension
        elif self.type == 'pmc':
            self.file_name = self.pmc_id + self.file_extension
            
        self.path = os.path.join(self.dump_dir, self.file_name)

        if not os.path.isdir(self.dump_dir):
            os.makedirs(self.dump_dir)


        (self.file, self.text, self.file_tree) = self._load_file(debug=debug)
        self.xml = et.tostring(self.file_tree)


    def _load_file(self, debug=False):
        '''Loads file object from storage. If it is not available attempts download from internet.'''
        #if not self.file_name in os.listdir(self.dump_dir):
        if not os.path.exists(self.path):
            self._download_xml()
        else: 
            if debug == True:
                print('File Available', self.file_name)

        file = codecs.open(self.path, 'rb', 'utf-8')
        #self.text = pickle.load(self.file)
        text = file.read()
        file.close()
        file_tree = et.fromstring(text)
        return (file, text, file_tree)

    def _download_xml(self):
        '''Download data from Pubmed Central (Entrez).'''

        if self.entrez_email == None:
            print('DOWNLOAD OF NEW DATA CANNOT BE COMPLETED, PLEASE RE-INSTANTIATE DUMPFILE OBJECT PROVIDING YOUR EMAIL.')
            return None

        Entrez.email = self.entrez_email
        print('ATTEMPT TO DOWNLOAD DATA FOR PubmedID', self.pmid, '--', 'PMCID ', self.pmc_id)

        if self.type == 'pmc' and self.pmc_id != None:
            handle = Entrez.efetch(db='pmc', id=self.pmc_id, retmode="xml")
        elif self.type == 'pubmed' and self.pmid != None:
            handle = Entrez.efetch(db="pubmed", id=self.pmid, retmode="xml")
        else:
            raise ValueError

        file_list = []
        tree = et.parse(handle)
        # hand generator to lxml
        
        download_file = open(self.path, 'w', encoding='utf-8')
        download_file.write(et.tostring(tree, encoding='unicode', pretty_print=True))
        download_file.close()

    def get_tree(self):
        return self.file_tree
        
    def write_tree(self, out_path):
        file_element_tree = et.ElementTree(self.file_tree)
        file_element_tree.write(out_path, pretty_print=True)
        print('XML file saved to', out_path)



class PMID_dump_file(DumpFile):

    def __init__(self, doc_id, dump_dir, id_mapper=None, entrez_email=None, debug=False):
        file_extension = '.pxml'
        super().__init__(doc_id, dump_dir, file_extension, id_mapper=id_mapper, entrez_email=entrez_email, debug=debug)

        self.section_list = ['title', 'abstract']
        self.title_text_dict = None

        self.subheadings = None 
            
    def clean_tree(self, tree='default'):
        '''
        Removes in-text tags that are not relevant and that cause the xml parser to break
        TO DO: adapt to PXML
        '''
        if tree == 'default':
            tree = copy.deepcopy(self.file_tree)  # use the whole file tree as a default value
        else:
            tree = copy.deepcopy(tree)

        if tree is None:
            print('ERROR, NO DATA FOR:', self.pmid)
            return None
        
        present_tags = list(set([element.tag for element in tree.iter()]))
        #print('Present Tags', present_tags)
        strip_tags = ['br', 'i', 'sup', 'b', 'sub']
        #print 'Tags to be stripped: ', ', '.join(strip_tags)
        # Remove in-text tags
        for strip_tag in strip_tags:
            et.strip_tags(tree, strip_tag)
        return tree

    def get_language(self):
        tree = self.file_tree
        #print(et.tostring(tree, pretty_print=True))

        language_element = tree.find('.//Language')
        if language_element is not None:
            return language_element.text
        else:
            return None

    def get_country(self):
        tree = self.file_tree
        #print(et.tostring(tree, pretty_print=True))

        country_element = tree.find('.//Country')
        if country_element is not None:
            return country_element.text
        else:
            return None

    def publication_types(self):
        tree = self.file_tree
        pub_types_element = tree.find('.//PublicationTypeList')
        pub_type_list = list()
        for pub_type_info in pub_types_element:  # PublicationType
            pub_type_list.append(pub_type_info.text)

        return tuple(pub_type_list)


    def get_date(self):
        date_field = self.file_tree.find('.//PubDate')
        if date_field is None:
            return None
        else:
            year = date_field.find('.//Year')
            if year is not None:
                return year.text

            else:  # use medline date
                year_re = re.compile('[0-9]{4}')
                medline_date = date_field.find('.//MedlineDate')
                if medline_date is None:
                    return None
                medline_date_text = medline_date.text
                year_list = year_re.findall(medline_date_text)
                return year_list[0]
            return None


    def get_title_text_dict(self, subheadings=True):
        '''
        Parse xml and return a dictionary containing the title and text only.
        subheadings are included by default into the text followed by a colon.
        '''

        if subheadings != self.subheadings:  # re-make dictionary with different subheading mode
            self.title_text_dict = None
            self.subheadings = subheadings

        if self.title_text_dict is not None:
            return self.title_text_dict

        title_text_dict = dict()
    
        tree = self.file_tree
        #print(et.tostring(tree, pretty_print=True))

        title_entry = tree.find('.//ArticleTitle')
        if title_entry is None:
            title_entry = tree.find('.//BookTitle')

        title = self.clean_tree(title_entry)
        if title is None:
            print('No title')
            return None
        title_text_dict['title'] = title.text

        abstract_list = []  # add all parts of an abstract
        for abs_entry in tree.findall('.//Abstract'):
            abs_tree = self.clean_tree(abs_entry)
            if abs_tree is None:  # no abstract found
                continue

            for abs_part in abs_tree.iterfind('.//AbstractText'):
                attribs = abs_part.attrib
                try:
                    sec_title = attribs['Label']
                except KeyError:
                    sec_title = None

                sec_text = abs_part.text

                if subheadings is True and sec_title is not None:
                    if sec_text:
                        # normalize unicode of section text 
                        # (especially important for normalization of a range of different possible whitespace characters)
                        sec_text = unicodedata.normalize('NFKD', sec_text)
                        abstract_list.append(sec_title + ': ' + sec_text)
                    else:
                        abstract_list.append(sec_title)
                else:
                    abstract_list.append(sec_text)
            
        abstract_string = ' '.join([a for a in abstract_list if a])
        # removing br-tags #is this still necessary?
        abstract_string = re.sub('(<br><br>|<br />)', '', abstract_string)
        title_text_dict['abstract'] = abstract_string

        self.title_text_dict = title_text_dict
        return title_text_dict


    def get_text_list(self, subheadings=True):
        ttdict = self.get_title_text_dict(subheadings=subheadings)
        sec_list = ['title', 'abstract']
        text_list = [ttdict[k] for k in sec_list]
        return text_list

    def get_abstract_sections(self):
        abstract_sec_dict = dict()

        for abs_part in self.file_tree.iterfind('.//AbstractText'):
                attribs = abs_part.attrib
                if 'Label' in attribs:
                    sec_title = attribs['Label']
                else:
                    sec_title = None
                if sec_title:
                    sec_text = abs_part.text
                    abstract_sec_dict[sec_title] = sec_text

        return abstract_sec_dict

    def has_abstract_text(self):
        '''
        Check if an document has an associated abstract text;
        Returns True only for documents that do not only have a title but also an abstract text
        '''
        title_text_dict = self.get_title_text_dict()
        if not title_text_dict:
            return False
        if title_text_dict['abstract']:
            return True
        else:
            return False

    def get_mesh(self):
        '''
        A list of MeSH headings/descriptors for the current record.
        Associated qualifiers/subheadings can be accessed with mesh_heading.qualifiers
        '''
        mesh_list = list()
        heading_section = self.file_tree.iterfind('.//MeshHeading')
        if heading_section == None:
            return mesh_list
        for head_item in heading_section:
            mesh_list.append(MeshHeading(head_item))
        return mesh_list

    def export_bioc(self, out_path, replace_abbrevs=False, include_offsets=True, subheadings=False, debug=False):
        from bioc_file import BioCFile, BioCCollection, BioCDocument, BioCPassage, BioCInfon
        from bioc_tools import OffsetCount
        

        ttc = self.get_title_text_dict(subheadings=subheadings)

        if replace_abbrevs is True:
            abbrev_matcher = AbbrevMatcher()
            for title, text_field in ttc.items():
                text_field_new = abbrev_matcher.replace_abbrevs(text_field)
                ttc[title] = text_field_new

        bioc_collection = BioCCollection()
        bioc_doc = BioCDocument(id=self.pmid)

        if include_offsets is True:
            offset_count = OffsetCount()
        else:
            offset_count = None

        for sec_name in self.section_list:
            if sec_name in ttc:
                section_text = ttc[sec_name]
                infons_dict = {'type': sec_name}
                bioc_passage = BioCPassage(infons_dict=infons_dict, content=None, text=section_text, offset=offset_count)
                # to do : include offsets!
                bioc_doc.add_passage(bioc_passage)

                if include_offsets is True:
                    offset_count.update(bioc_passage.length)

        bioc_collection.add_document(bioc_doc)

        file_out = BioCFile(bioc_collection=bioc_collection)
        file_out.write_bioc(out_path)


        
        
class MeshHeading(object):
    '''
    One MeSH descriptor/heading.
    self.qualifiers: list of associated qualifiers
    self.list: list of descriptor/qualifier pairs
    '''
    def __init__(self, mesh_element):
        self.alllist = list()
        for item in mesh_element:
            self.alllist.append(MeshDescQual(item))

        descriptors = [mh for mh in self.alllist if mh.type == 'Descriptor']
        if len(descriptors) > 1:
            print('Error: more than one MeSH descriptor found in Heading Definition')
            raise Exception
        else:
            self.descriptor = descriptors[0]

        self.qualifiers = [mh for mh in self.alllist if mh.type == 'Qualifier']

        self.list = [str(self.descriptor) + '/' + str(qual.str) for qual in self.qualifiers]
        self.ui = self.descriptor.ui
        self.type = self.descriptor.type
        self.str = self.descriptor.str

    def __repr__(self):
        return self.descriptor

    def __iter__(self):
        for qual in self.qualifiers:
            yield qual

    def has_qualifiers(self):
        if self.qualifiers != []:
            return True
        return False


class MeshDescQual(object):
    '''
    One MeSH descriptor/heading or qualifier/subheading.
    type: Descriptor or Qualifier
    ui: Unique Identifier (MeSH ID)
    majr: Major Topic of Abstract (Y/N)
    '''
    def __init__(self, heading_element):
        self.tag = heading_element.tag
        self.type = self.tag.rstrip('Name')
        self.ui = heading_element.get('UI')  # Unique Identifier (MeSH ID)
        self.majr = heading_element.get('MajorTopicYN')
        self.str = heading_element.text

    def __repr__(self):
        return self.str

    def __str__(self):
        return self.str


class TitleContentDict(object):
    def __init__(self, tree, include_titles=False, include_references=False, include_additional_abstracts=False, debug=False):
        self.tree = tree
        self.data = {}
        self.current_section = None
        self.section_tags = ['article-title', 'abstract', 'body']
        self.title_tags = ['title']
        # caption: avoids inclusion of irrelevant tags (and tag content)
        self.continue_tags = ['label', 'caption', 'dips-formula', 'inline-formula', 'xref', 'author-notes', 'related-article']
        self.titles = []    # section titles in the right order
        self.include_references = include_references
        self.include_titles = include_titles
        self.include_additional_abstracts = include_additional_abstracts #include abstracts that go beyond the main abstract
        self.debug = debug


    def process_text_element(self, text_element, remove_sb=True, replace_newline=True, replace_multiple_ws=True):
        for sel in self.continue_tags:
            et.strip_tags(text_element, sel)

        text = text_element.text

        processed_text = self.process_text(
            text,
            remove_sb=remove_sb,
            replace_newline=replace_newline,
            replace_multiple_ws=replace_multiple_ws
                                    )

        return processed_text
        
    def process_text(self, text, remove_sb=True, replace_newline=True, replace_multiple_ws=True):
        '''
        Args:
            remove_sb: 
                remove empty or almost empty (containing only certain punctuation) squared brackets 
                (originating from references) from text.
            replace_newline: 
                replace all newline characters by single whitespace characters
            replace_multiple_ws: 
                replace all occurrences of 2 or more sequential whitespace characters by 
                single whitespace characters
        '''

        if text is None:
            return ''
        if replace_newline is True:
            text = text.replace('\n', ' ')
        if remove_sb is True:
            rep_re = '\[[\s;,\./-]*\]'
            text = re.sub(rep_re, '', text)
        if replace_multiple_ws is True:
            text = re.sub('\s{2,}', ' ', text)
        text = text.strip()
        text = unicodedata.normalize('NFKD', text)
        text = text.strip()
        #print([text])
        return text
        
    def parse_data(self, tree='default', current_title=''):
        seen_sections = list()  # keep record of seen sections

        if tree == 'default':
            tree = self.tree

        for child in tree.getchildren():  # only get immediate children
            if child.tag in self.continue_tags:
                # print('continue tag:', child.tag)
                continue  # ignore current element and all its children elements
            if self.include_references is False and child.tag in self.titles:
                # the child's tag is already in the list of processed section headings
                if child.tag == 'article-title':
                    if 'body' in seen_sections:
                        break  # it is most likely a title of a reference (which are not to be included)
                    else:
                        #print(child.text)
                        # e.g. abstracts of related papers
                        pass

            if child.tag in self.section_tags and child.tag not in seen_sections:  # a new section starts
                #print(child.tag)
                #print('abstract type', child.get("abstract-type"))
                seen_sections.append(child.tag)
                self.current_section = child.tag
                current_title = self.process_text(self.current_section)
                # add the mapping between the current title and the text to the data
                self.data[current_title] = [self.process_text_element(child)]
                self.titles.append(current_title)
    

            elif child.tag in self.title_tags:
                current_title = current_title + ':' + self.process_text_element(child)
                
                if self.include_titles is True:  # insert the title as part of the text
                    self.data[current_title] = [self.process_text_element(child) + ': ']
                else:
                    self.data[current_title] = []
                    
                self.titles.append(current_title)

                if self.include_references is False and child.text == 'References':
                    break  # everything after 'References' is not to be included

            elif child.tag == 'article-meta' and 'abstract' not in seen_sections:
                # in some cases 'abstract' might be in a 'article-meta' tag and this tag should be parsed
                pass

            else:
                try:
                    self.data[current_title].append(self.process_text_element(child))
                except KeyError:
                    if self.debug is True:
                        print('Warning: Unbound content found')  # content cannot be associated to any part of the article
                        print(child.text)

                    self.data[current_title] = [self.process_text_element(child)]  # in case the title is (still) empty

            self.parse_data(tree=child, current_title=current_title)  # recursively parse everything
            
    def iter_content(self):
        data_dict = self.get_dict()
        for title in self.titles:
            yield (title, data_dict[title])
            
    def get_dict(self):
        if not self.data:
            self.parse_data()
        tcd = dict()
        for title, content_list in self.data.items():
            tcd[title] = ' '.join(content_list).strip()
        return tcd

    def get_title_text_dict(self, include_titles=False, include_newlines=False):
        '''
        Generate dictionary with the following keys corresponding to the core titles of a fulltext paper:
        'article-tilte', 'abstract', 'body'
        '''
        core_titles = ['article-title', 'abstract', 'body']
        if not self.data:
            self.parse_data()
        data_dict = self.get_dict()
        #print('DATA DICT', data_dict)
        ttd = dict()
        for title in self.titles:
            content = data_dict[title]

            title_list = title.split(':')
            current_title = title_list[0]

            if include_newlines is True:  # include newlines after every title or content section
                title = title + ': \n'
                content = content + '\n'

            if current_title not in ttd:
                ttd[current_title] = [content]
            else:
                if include_titles:
                    ttd[current_title].append(title)
                ttd[current_title].append(content)

        for ct in core_titles:  # add core titles that are not yet present (e.g. if an article is missing an abstract)

            if ct not in ttd:
                if self.debug:
                    print('CORE TITLE NOT FOUND', ct)
                ttd[ct] = ''
            else:
                ttd[ct] = ' '.join(ttd[ct]).strip()
        return ttd

    def get_abstract_types(self):
        return [i.get('abstract-type') for i in self.tree.iterfind('abstract')]


class PMC_dump_file(DumpFile):

    def __init__(self, doc_id, dump_dir, id_mapper=None, entrez_email=None, debug=False):
        file_extension = '.nxml'
        super().__init__(doc_id, dump_dir, file_extension, id_mapper=id_mapper, entrez_email=entrez_email, debug=debug)
        
        if self.debug:
            print('MAPPED IDs: ', self.pmid, '--', self.pmc_id)


    def fulltext_available(self):
        '''
        Checks if fulltext is available (evidence through comment in xml) and return True/False.
        Comment in xml: <!--The publisher of this article does not allow downloading of the full text in XML form.-->
        '''
        no_fulltext_comment = '<!--The publisher of this article does not allow downloading of the full text in XML form.-->'
        comment_list = [str(c) for c in self.file_tree.iter(tag=et.Comment)]
        if no_fulltext_comment in comment_list:
            return False
        else:
            return True
       
    def clean_tree(self, tree='default', removables=['license', 'notes', 'trans-abstract', 'trans-title']):
        '''
        Removes in-text tags that are not relevant and that cause the xml parser to break

        Args:
            tree: xml tree; default is the instance's whole file tree
            removables: list of tags to be discarded together with their contents
        '''
        if tree == 'default':
            #tree = self.file_tree
            tree = copy.deepcopy(self.file_tree)
        else:
            tree = copy.deepcopy(tree)

        if tree is None:
            print('ERROR, NO DATA FOR:', self.pmc_id)
            return None
        
        core_relevant_tags = ['body', 'p', 'sec', 'label', 'title', 'abstract', 'article-title']
        # some of them need to be present in order to be able to remove them later on
        additional_relevant_tags = [
            'list',
            'boxed-text',
            'list-item',
            'caption',
            'inline-formula',
            'author-notes',
            'article-meta'
                                    ]

        relevant_tags = core_relevant_tags + additional_relevant_tags
        present_tags = list(set([element.tag for element in tree.iter()]))

        # define tags to remove together with all children
        deltag = 'remove-tag'
        remove_tags = [deltag] + removables
        # deltag is a dummy tag for which is used to replace all tags with specific attribs to be removed
        remove_tags_with_attribs = [
            ('xref', 'ref-type="bibr"'),  # bibligraphic references (mostly numbers)
            ('xref', 'ref-type="contrib"'),  # author contributions
            ('xref', 'ref-type="aff"'),  # affiliations
            ('xref', 'ref-type="fn"'),  # reference to footnote (mostly numbers)
            ('xref', 'ref-type="bio"'),
            ('xref', 'ref-type="author-notes"')
                   
                                    ]  
        # xml tags to be stripped from the input xml
        strip_tags = [tag for tag in present_tags if tag not in relevant_tags]  
        if self.debug:
            print('STRIP TAGS', strip_tags)
            print('REMOVE TAGS', remove_tags)
            print('REMOVE TAGS WITH ATTRIBS', remove_tags_with_attribs)

        # for el in tree.iterfind('.//xref'):
        #     print(el.get('ref-type'))

        # prepare specified tags with a specific attribute to be removed by converting them to an identifiable dummy tag
        tags_attribs_findstrings = [".//{tag}[@{attrib}]".format(tag=tag, attrib=attrib)
                                    for (tag, attrib) in remove_tags_with_attribs]

        for fs in tags_attribs_findstrings:
            for el in tree.iterfind(fs):
                el.tag = deltag

        for remove_tag in remove_tags:
            et.strip_elements(tree, remove_tag, with_tail=False)

        # Remove in-text tags
        for strip_tag in strip_tags:
            et.strip_tags(tree, strip_tag)
            
        return tree


    def get_abstract_xml(self, remove_newline=True):
        '''Find abstract section and return a dictionary mapping abstract section headings to their content.'''
        tree = self.clean_tree(self.file_tree)
        if tree is None:
            print('ERROR, No Data')
            return dict()
        abstract_tree = [et.tostring(i, pretty_print=True, encoding='unicode') for i in tree.iterfind(".//abstract")]
        return abstract_tree

    def get_title(self):
        '''Get text version of the title of the article.'''
        article_title = self.file_tree.find(".//article-title")
        if article_title is not None:
            return article_title.text
        else:
            return None

    @staticmethod
    def get_sec_title(sec, remove_newline=True):
            sec_title = sec.find(".//title")
            if sec_title is None:
                return None
            sec_title = sec_title.text
            if remove_newline is True:
                sec_title = sec_title.replace("\n", "")
                sec_title = unicodedata.normalize('NFKD', sec_title)
            return sec_title

    def get_text(sec, replace_newline=True, replace_multiple_ws=True):
        '''
        Args:
            replace_newline: 
                replace all newline characters by single whitespace characters
            replace_multiple_ws:
                replace all occurrences of 2 or more sequential whitespace characters by
                single whitespace characters
        '''
        if sec is None:
            return ''
        if sec.text:
            sec_text = sec.text
            if replace_newline is True:
                sec_text = sec_text.replace('\n', '')
            if replace_multiple_ws is True:
                sec_text = re.sub('\s{2,}', ' ', sec_text)
            return sec_text
        else:
            return ''

    @staticmethod
    def get_sec_content(sec, replace_newline=True, replace_multiple_ws=True):
        '''
        Args:
            replace_newline:
                replace all newline characters by single whitespace characters (default)
            replace_multiple_ws:
                replace all occurrences of 2 or more sequential whitespace characters
                by single whitespace characters (default)
        '''
        if sec is None:
            return ''
        sec_children = sec.getchildren()
        for child in sec_children:
            if child.tag == "p":
                content = child.text
                if replace_newline is True:
                    content = content.replace("\n", " ")
                if replace_multiple_ws is True:
                    content = re.sub('\s{2,}', ' ', content)
                content = unicodedata.normalize('NFKD', content)
                return content
            elif child.tag == "sec":  # no content before next section
                return ''

    def transform_to_list(self, tree='default', include_references=False):
        '''Takes an xml object and transforms all its content to a list of text.'''
        clean_tree = self.clean_tree(tree=tree)
        text_list = list()
        continue_tags = ['pmc-articleset', 'label']  # nothing happens for these tags
        for i in clean_tree.iter():
            if i.tag in continue_tags:
                continue
            if i.text == 'References' and include_references is False:
                break  # do not include reference section or anything following it
            elif i.tag == 'title':
                text = i.text + ': '
            else:
                text = i.text
            if text:
                text = text.replace('\n', '')
                text_list.append(unicodedata.normalize('NFKD', text))
        return text_list


    @staticmethod
    def get_subsection_headings(sec):
        subsections = [i for i in sec if i.tag == 'sec']
        return [PMC_dump_file.get_sec_title(i) for i in subsections]

    def body_section_headings(self):
        ''' Returns all first level section headings '''
        article_body = self.file_tree.find(".//body")
        if article_body is None:
            return None  # no body available
        return PMC_dump_file.get_subsection_headings(article_body)  # first level section heading


    def get_date(self):
        date_field = self.file_tree.find('.//pub-date[@pub-type="epub"]')
        if date_field is None:
            return None
        else:
            year = date_field.find('.//year')
            return year.text


    def title_content_dict(self, sec, remove_newline=True):
        '''
        Goes through the xml and generates a dictionary structure mapping section headings
        to section content in a recursive way.
        Dictionary keys on each level:
        'abstract' (only first level), 'body' (only first level), 'text', 'subsections';
        'subsections' contains a dictionary
        * TODO: PROBLEM with sections without a title! Sections without a title are lost. [CHANGE THIS!!]
        '''
        # TODO: replace this using TitleContentDict object
        raise NotImplementedError

    def get_title_text_dict(self, tree='default', include_titles=False, include_newlines=False, include_references=False):
        if tree == 'default':
            clean_tree = self.clean_tree()
        else:
            clean_tree = self.clean_tree(tree=tree)
        title_content = TitleContentDict(clean_tree, include_references=include_references)
        return title_content.get_title_text_dict(include_titles=include_titles, include_newlines=include_newlines)

    def get_title_content(self, include_references=False):
        clean_tree = self.clean_tree()
        title_content = TitleContentDict(clean_tree, include_references=include_references)
        return title_content

    def get_title_sections_tuples(self):
        clean_tree = self.clean_tree()
        title_content = TitleContentDict(clean_tree)
        return title_content.iter_content()

    def get_text_list(self, subheadings=False):
        ttdict = self.get_title_text_dict()
        sec_list = ['article-title', 'abstract', 'body']
        text_list = [ttdict[k] for k in sec_list]
        return text_list


    def export_bioc(self, out_path, replace_abbrevs=False, include_offsets=True, debug=False):
        from bioc_file import BioCFile, BioCCollection, BioCDocument, BioCPassage, BioCInfon
        from bioc_tools import OffsetCount

        clean_tree = self.clean_tree()

        bioc_collection = BioCCollection()
        bioc_doc = BioCDocument(id=self.doc_id)

        if include_offsets is True:
            offset_count = OffsetCount()
        else:
            offset_count = None


        title_content = TitleContentDict(clean_tree)

        if replace_abbrevs is True:
            ttc_new = dict()
            abbrev_matcher = AbbrevMatcher(debug=debug)
            for title, textfield in title_content.iter_content():
                new_text = abbrev_matcher.replace_abbrevs(textfield)
                ttc_new[title] = new_text


        for title, text in title_content.iter_content():
            if replace_abbrevs is True:
                text = ttc_new[title]

            infons_dict = {'type': title}
            bioc_passage = BioCPassage(infons_dict=infons_dict, content=None, text=text, offset=offset_count)
            # to do : include offsets!
            bioc_doc.add_passage(bioc_passage)

            if include_offsets is True:
                offset_count.update(bioc_passage.length)

        bioc_collection.add_document(bioc_doc)

        file_out = BioCFile(bioc_collection=bioc_collection)
        file_out.write_bioc(out_path)

    def check_for_body(self):
        ttd = self.get_title_text_dict()
        if ttd['body']:
            return True
        else:
            return False

    def get_xml_lang(self):
        '''
        Get xml language from xml:lang attribute
        *TODO: problem - this does not work anymore after fulltext detection; does self.file_tree change?
        '''
        tree = self.file_tree
        article = tree.find(".//article")
        if article is None:
            print('No article for', self.pmc_id)
            print(et.tostring(tree))
            return 'no_article'

        try:
            lang = article.attrib['{http://www.w3.org/XML/1998/namespace}lang']
            return lang
        except KeyError:
            return 'no_lang_attrib'


    def get_abstract_types(self):
        tree = self.file_tree
        return [i.get('abstract-type') for i in tree.iterfind('.//abstract')]
        


                

def main():
    """
    Invoke this module as a script
    """
    description="Module for downloading, storing, loading and handling of Pubmed and PMC dump files in XML format."
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('-l', '--logfile', dest='logfilename',
                        help='write log to FILE', metavar='FILE')
    parser.add_argument('-q', '--quiet',
                        action='store_true', dest='quiet', default=False,
                        help='do not print status messages to stderr')
    parser.add_argument('-d', '--debug',
                        action='store_true', dest='debug', default=False,
                        help='print debug information')
    
    
    
    args = parser.parse_args()
    print('Args:', args)
    
    pass
    
if __name__ == '__main__':
    main()