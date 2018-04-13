#!/usr/bin/env python

"""
Module for finding abbreviations in text (if they have the standard format of a full form followed by 
the abbreviation in brackets).
Contains methods for extracting abbreviations from a file or from a collection and replacing abbreviations
in text by full forms.

Author: Tilia Ellendorff

"""

import sys
import os
import argparse
import string
import re


sys.path.append('lib')


class AbbrevMatcher(object):
    def __init__(self, min_len=2, max_len=4, long_version_max=7, abbrev_exclude=list(), debug=False):
        self.data = dict()
        self.bracket_match = re.compile('\s?\([^\)^\(\)]+\)[^\w]')
        self.min_len = min_len
        self.max_len = max_len

        self.long_version_max = long_version_max

        self.abbrev_exclude = abbrev_exclude

        self.debug = debug

    def prefilter_candidate(self, abbrev_candidate):
        '''
        Prefilter abbreviation candidates based on features related to their own string.
        candidates are sorted out:
        - if their length is longer than the maximum length
        - if they contain a '='
        - if they only contain non-alphabetical characters
        '''

        # rule: needs to be of defined length
        if not self.min_len <= len(abbrev_candidate) <= self.max_len:
            return False
        if '=' in abbrev_candidate:  # rule: no '=' allowed
            return False
        # rule: needs to contain at least one alphabetic character
        if not any(char.isalpha() for char in abbrev_candidate):
            return False
        # rule: exclude all lower-case abbreviations
        if abbrev_candidate.islower():
            return False
        if abbrev_candidate in self.abbrev_exclude:
            return False
        return True

    def prefilter_long_version(self, long_version_candidate, abbrev_candidate):
        if len(long_version_candidate) > 50:
            return False
        if ' .' in long_version_candidate:
            return False
        num_punct_chars = len([c for c in long_version_candidate if c in string.punctuation])
        if num_punct_chars > 2:
            return False
        # rule: abbreviation should not be part of the long version
        if abbrev_candidate.lower() in long_version_candidate.lower():
            return False
        # rule: all characters of the abbreviation candidate should be part of the long version
        if not all(char.lower() in long_version_candidate.lower() for char in abbrev_candidate):
            return False
        return True


    def filter_long_version_list(self, long_version_list, abbrev_candidate):
        lv_list = list()
        if len(set([lv.lower() for lv in long_version_list])) == 1:
            return [long_version_list[0]]
        if len(set([lv.lower().strip('gene') for lv in long_version_list])) == 1:
            lv_list.extend(lv for lv in long_version_list if not lv.endswith('gene'))
            return lv_list
        for lv in long_version_list:
            lv_starts = [lt[0].lower() for lt in lv.split()]
            if all(char.lower() in lv_starts for char in abbrev_candidate):
                lv_list.append(lv)
                return lv_list

            elif not set(string.punctuation).intersection(lv):
                lv_list.append(lv)

        if lv_list:
            return lv_list
        else:
            return long_version_list

    def match_abbrevs(self, text_string):
        '''Uses regular expressions for matching abbreviations in the text'''

        # bracket_match = '\([^\)]+\)'
        # m = re.findall(bracket_match, text_string)
        # if m:
        #   print(text_string)
        #   print(m)

        stopwords = ['and', 'also']

        
        for m in self.bracket_match.finditer(text_string):
            abbrev_candidate = m.group().strip().rstrip(string.punctuation).rstrip(')').lstrip('(')
            if self.prefilter_candidate(abbrev_candidate) is False:
                continue
            else:
                #print(m.start(), m.group())

                start_char = abbrev_candidate[0]
                abbrev_offset = m.start()

                # only look at text string before abbreviation
                pre_string = text_string[:abbrev_offset]  
                pre_string_list = re.split('[\/\s-]', pre_string)
                # reverse list of tokens before abbreviation (start right before abbreviation candidate)
                pre_string_list_rev = [w for w in reversed(pre_string_list)]
                #print(pre_string_list_rev)
                
                ex_list = list()  # tokens of possible long version candidate

                for pos, w in enumerate(pre_string_list_rev):  # pos is the position before the abbreviation candidate
                    
                    if pos == 0:  # skipping first word directly before
                        try:
                            first_before_start = w[0].lower()  # token immediately before abbreviation
                            first_before = w
                        except IndexError:
                            #print(w, 'no chars')
                            #print(text_string)
                            first_before_start = None
                        ex_list.append(w)
                        continue

                    # strip first token of long version candidate from punctuation characters
                    w_test = w.strip(string.punctuation)
                    try:
                        word_start = w_test[0]  # first character of the word
                    except IndexError:
                        continue
                    if word_start.lower() == start_char.lower():
                        ex_list.append(w)
                        if not w.lower() in stopwords:
                            break  # first token of long version found

                    # it is 8 tokens long already and the token immediately before has the same start character
                    if pos >= 7 and first_before_start == start_char.lower():
                        ex_list = [pre_string_list_rev[0]]
                        break  # use token immeadiately before abbreviation as first token in long version
                    if pos == len(pre_string_list_rev) - 1:  # if end of prestringlist is reached
                        if pos <= self.long_version_max:
                            ex_list = [pre_string_list_rev[0]]
                        else:
                            ex_list = list()  # delete everything if string gets too long
                    else:
                        ex_list.append(w)
                
                if not ex_list == []:
                    long_version = ' '.join([w for w in reversed(ex_list)]).strip(string.punctuation)

                    # long version has been filtered based on its string features
                    if self.prefilter_long_version(long_version, abbrev_candidate) is False:
                        continue

                    if abbrev_candidate not in self.data:
                        self.data[abbrev_candidate] = [long_version]
                    else:
                        self.data[abbrev_candidate].append(long_version)
                    #print(long_version)

        return self.data

    def clear_cache(self):
        self.data = dict()

    def replace_abbrevs(self, instring):
        if instring is None:
            return None
            
        self.match_abbrevs(instring)
        for abbrev, long_versions in self.data.items():
            if len(set([lv.lower() for lv in long_versions])) == 1:
                instring = instring.replace(abbrev, long_versions[0])
            else:
                filtered_long_versions = self.filter_long_version_list(long_versions, abbrev)
                if len(filtered_long_versions) == 1:
                    instring = instring.replace(abbrev, filtered_long_versions[0])
                else:
                    instring = instring.replace(abbrev, long_versions[0])
                    if self.debug is True:
                        print('Ambiguous abbreviation')
                        print(abbrev, long_versions)
                        #raise ValueError
        return instring

    def extract_from_file(self, file_path):
        '''
        Finds abbreviations in one text file.
        And returns a dictionary mapping abbreviations to their long forms.
        '''
        with open(file_path, 'r') as f:
            f_lines = f.readlines()
            for line in f_lines:
                self.match_abbrevs(line)
        return self.data


    def extract_from_collection_files(self, in_dir, doc_ids=list()):
        if doc_ids == []:
            filenames = [f for f in os.listdir(in_dir) if f.endswith('.txt')]
            print(len(filenames), 'files will be considered')
            for fn in filenames:
                file_path = os.path.join(in_dir, fn)
                self.process_file(file_path)  # update the collection abbreviation dictionary
        else:
            for doc_id in doc_ids:
                file_path = os.path.join(in_dir, doc_id + '.txt')
                self.process_file(file_path)

        # print('Discovered Abbreviations')
        # for abbrev, ex_list in collection_data.items():
        #   print(abbrev, ':', list(set(ex_list)))
        print(len(self.data), 'Abbreviations Discovered')
        return self.data


def process(args=None):
    pass
    
    
    
def main():
    """
    Invoke this module as a script
    """
    description = "...."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-l', '--logfile', dest='logfilename',
                        help='write log to FILE', metavar='FILE')
    parser.add_argument('-q', '--quiet',
                        action='store_true', dest='quiet', default=False,
                        help='do not print status messages to stderr')
    parser.add_argument('-d', '--debug',
                        action='store_true', dest='debug', default=False,
                        help='print debug information')
    
    parser.add_argument('--in_path', metavar='in_path', type=str, required=True,
                        help='Path to bioc directory/file to be converted.')
    parser.add_argument('--out_path', metavar='out_path', type=str, required=True,
                        help='Path to output file/directory.')
    
    
    
    args = parser.parse_args()
    print('Args:', args)
    
    process(args=args)

    sys.exit(0)  # Everything went ok!
    
    
if __name__ == '__main__':
    main()
