# dictionary.py
# Python version of WatDiv/dictionary
# T. Masuda, 2023/12/1

from enum import Enum, auto
# #include "dictionary.h"
#
# #include <algorithm>
# #include <fstream>
# #include <iostream>
# #include <sstream>

# using namespace std;

# namespace DICTIONARY_TYPES {
#     enum enum_t {
#         ENGLISH_WORDS,
#         FIRST_NAMES,
#         LAST_NAMES,
#         COUNT
#     };
# };
class DICTIONARY_TYPES(Enum):
    ENGLISH_WORDS = auto
    FIRST_NAMES = auto
    LAST_NAMES = auto


DICTIONARY_TYPES_COUNT = len(DICTIONARY_TYPES)


class Dictionary:
    """

    Attributes:

    """
#         vector<string>          _words[(int)DICTIONARY_TYPES::COUNT];
    _words: list[list[str]] = [[], [], []]
# dictionary * dictionary::_instance = NULL;
    _instance = None

# dictionary::dictionary(){
# }

# dictionary::~dictionary(){
# }

# dictionary * dictionary::get_instance(){
    @staticmethod
    def get_instance():
        """

        Args:

        Returns:

        """
#     if (_instance==NULL){
        if Dictionary._instance == None:
#         _instance = new dictionary();
            Dictionary._instance = Dictionary()
#     }
#     return _instance;
        return Dictionary._instance
# }
        pass

# void dictionary::destroy_instance(){
#     if (_instance!=NULL){
#         delete _instance;
#     }
#     _instance = NULL;
# }

# void dictionary::init (const char * words_filename, const char * firstnames_filename, const char * lastnames_filename){
    def init(self, words_filename: str, firstnames_filename: str, lastnames_filename: str):
        """

        Args:

        Returns:

        """
#     //cout<<"Initializing dictionary..."<<"\n";
        print("Initializing dictionary...")

#     for (int i=0; i<DICTIONARY_TYPES::COUNT; i++){
        for i in range(DICTIONARY_TYPES_COUNT):
#         ifstream fis;
#         switch (i){
#             case 0:{
            if i == 0:
#                 fis.open(words_filename);
                file_open = words_filename
#                 break;
#             }
#             case 1:{
            elif i == 1:
#                 fis.open(firstnames_filename);
                file_open = firstnames_filename
#                 break;
#             }
#             case 2:{
            else:
#                 fis.open(lastnames_filename);
                file_open = lastnames_filename
#                 break;
#             }
#         }
#         string line, token;
#         while (fis.good() && !fis.eof()){
            with open(file_open, 'r') as fis:
#             getline(fis, line);
                lines = fis.readlines()
                for line in lines:
    #             stringstream parser(line);
                    parser = line.replace('\n', '').split(' ')
    #             while (parser>>token){
                    for token in parser:
    #                 for (int k=0; k<token.size(); k++){
                        for k in range(len(token)):
    #                     if (((int)(token[k]))<32){
                            if ord(token[k]) < 32:
    #                         token[k] = (char) 32;
                                token = token[:k]+chr(32)+token[k+1:]
    #                     }
    #                     if (((int)(token[k]))>127){
                            if ord(token[k]) > 127:
    #                         token[k] = (char) 127;
                                token = token[:k]+chr(127)+token[k+1:]
    #                     }
    #                 }
    #                 _words[i].push_back(token);
                        Dictionary._words[i].append(token)
    #             }
#         }
                pass
#         fis.close();
#         sort(_words[i].begin(), _words[i].end());
            Dictionary._words[i].sort()
#     }
            pass
#     //cout<<"Dictionary initialized..."<<"\n";
        print("Dictionary initialized...")
# }
        pass  # end of dictionary::init

# unsigned int dictionary::word_count (DICTIONARY_TYPES::enum_t dictionary_type){
    def word_count(self, dictionary_type):
        """

        Args:

        Returns:

        """
#     return _words[((int)dictionary_type)].size();
        return len(self._words[int(dictionary_type)])
# }
        pass  # end of dictionary::word_count

# string * dictionary::get_word (DICTIONARY_TYPES::enum_t dictionary_type, unsigned int index){
    def get_word(self, dictionary_type, index):
        """

        Args:

        Returns:

        """
#     if (index<_words[((int)dictionary_type)].size()){
        if index < len(self._words[int(dictionary_type)]):
#         return &(_words[((int)dictionary_type)][index]);
            return self._words[int(dictionary_type)][index]
#     } else {
        else:
#         return NULL;
            return None
#     }
        pass  # end of if
# }
        pass  # end of dictionary::get_word

# pair<unsigned int, unsigned int> dictionary::get_interval (DICTIONARY_TYPES::enum_t dictionary_type, const string & range_min, const string & range_max){
    def get_interval(self, dictionary_type, range_min, range_max):
        """

        Args:

        Returns:

        """
#     // The interval is closed on the lower-bound and open on the upper-bound...
#     vector<string>::iterator low, high;
#     low = lower_bound(_words[((int)dictionary_type)].begin(), _words[((int)dictionary_type)].end(), range_min);
#     high = upper_bound(_words[((int)dictionary_type)].begin(), _words[((int)dictionary_type)].end(), range_max);
#     return pair<unsigned int, unsigned int> (low-_words[((int)dictionary_type)].begin(), high-_words[((int)dictionary_type)].begin());
        return (0, 0)
# }
