"""
version 1.0
author: Felix Plasser
description: Procedures for reading and writing files.
"""

import os, locale

class dict_plus(dict):
    """
    Extension of a dictionary where data can be read from an input file.
    """
    def __init__(self, file_name='', revert=False):
        if not file_name == '':
            if not revert:
                self.read_from_file(file_name)
            else:
                self.revert_read_from_file(file_name)
    
    def read_from_file(self,file_name):
        """
        Read data from a file for the dictionary.
        """
        r_file = open(file_name, 'r')
        for line in r_file:
            spaceind = line.find(' ')
            self[line[:spaceind]]=line[spaceind+1:-1]
            
    def revert_read_from_file(self,file_name):
        """
        Read data from a file for the dictionary. Keys and values are switched
        """
        r_file = open(file_name, 'r')
        for line in r_file:
            spaceind = line.find(' ')
            self[line[spaceind+1:-1]]=line[:spaceind]

def chmkdirs(dir):
    """
    Tries to change to the directory *dir* and creates it if it does not exist.
    """
    try:
        os.chdir(dir)
    except:
        os.makedirs(dir)
        os.chdir(dir)
        
def line_to_words(line, sep=[' '], convert=None):
    """
    For parsing files.
    Returns a list of words that were given in <line>, separated by <sep>.
    If line ends in'\n', this is rejected.
    If <not convert==None> values are converted according to convert
        e.g. <convert=float> for conversion into floating point numbers.
    -> it is probably more useful to use the string.split() function
    """
    if line == '':
        return []
    if line[-1] == '\n':
        line = line[:-1]
    ret_list = []
    tmp_str = ''
    for let in line:
        if not let in sep:
            tmp_str += let
        else:
            if not tmp_str == '':
                if not convert==None:
                    tmp_str = eval(convert+'(tmp_str)')
                ret_list.append(tmp_str)
            tmp_str = ''
    if not tmp_str == '':
        if not convert==None:
            tmp_str = eval(convert+'(tmp_str)')
        ret_list.append(tmp_str)
            
    return ret_list

def change_file_lines(file_name, ind_cont=[]):
    """
    Change lines in a file (lines are overwritten).
    <ind_cont> is a list with indices and the contents for this line, e.g. ind_cont=[[3,'test']] writes 'test' into line 3.
    """
    lines = open(file_name, 'r').readlines()
    for ind, cont in ind_cont:
        lines[ind-1]=cont + '\n'

    w_file = open(file_name, 'w')
    w_file.writelines(lines)
    w_file.close()
    
    
def get_file_keys(file_name, equ='=', sep=[' ',',']):
    """
    Parses statements in input files. All information is returned in a dictionary.
    "opt1=5" -> {"opt1":"5"}
    "opt2=1,2" -> {"opt2":["1","2"]}
    if <sep> is set to zero, then the output works according to original python syntax
    """
    out_dict = {}
    fil = open(file_name,'r')
    for line in fil:
        if line.count(equ) == 1:
            tmp = line_to_words(line,sep=sep+[equ])
            out_dict[tmp[0]]=tmp[1:]
        elif line.count(equ) > 1:
            print("line could not be parsed because there are more than one '%s'"%equ)
            print(line)
    
    fil.close()
    return out_dict
        
def to(start,end):
    """
    Interface to "range" which includes <end>.
    This is used to allow for more intuitive user input.
    """
    return range(start, end+1)

def list_print(list,format,ncols=8):
    """
    Print a list formatted.
    """
    outstr=''
    i=0
    for item in list:
        i+=1
        outstr+=format%item

        if i == ncols:
            print(outstr)
            i = 0
            outstr = ''

    print(outstr)

class table_maker:
    """
    Class for writing output in columns.
    It is probably more useful to use intrinsic python formatting instead.
    """
    def __init__(self, col_widths, cut=True, replace_list=[], digits=3):
        """
        Enter the widths of the columns in list <col_widths>.
        <replace_list> contains a double, list of items to be replaced,
            e.g. replace_list=[['.',',']] for European decimal notation.
        If <cut==True> content is cut to fit into the columns.
        <digits> optionally gives the number of digits for rounding.
        """
        self.col_widths = col_widths
        self.col_nr = len(self.col_widths)
        self.ret_string = '' # string to be returned
        self.cut = cut
        self.replace_list = replace_list
        self.digits = digits

    def write_line(self, words):
        """
        Writes a line with list <words> to the table.
        """    
        self.ret_string += self.ret_line(words) + '\n'

    def write_lines(self, words):
        """
        Writes lines with the appropriate number of columns to the table with a longer list words.
        """
        nr_lines = len(words) / self.col_nr
        if len(words) % self.col_nr != 0:
            nr_lines += 1        

        for i in xrange(nr_lines):
            #print words[i*self.col_nr : (i+1)*self.col_nr]
            self.write_line(words[i*self.col_nr : (i+1)*self.col_nr])
        
    def ret_line(self, words):
        """
        Return a formatted line with list <words>.
        """
        plus_string = ''
        for i, word in enumerate(words):
            
            if self.digits != None and type(word)==float:
                word = locale.format("% .*f", (self.digits, word))
            elif self.cut:
                word = str(word)[:(self.col_widths[i]-1)]
            
            plus_string += str(word).ljust(self.col_widths[i])

        for old,new in self.replace_list:
            plus_string = plus_string.replace(old, new)
    
        return plus_string
        
    def print_line(self,words):
        """
        Print the formatted line.
        """
        print(self.ret_line(words))
        
    def make_line(self, words):
        """
        Make a line with 
        """

    def return_table(self):
        """
        Return the table that has been written with write_line.
        """
        return self.ret_string

    def write_to_file(self, file_name):
        """
        Write the table to a file.
        """
        write_to_file(self.return_table(), file_name)
        
class csv_maker:
    """
    Class for making a csv file.
    """
    def __init__(self, sep=',', replace_list=[]):
        """
        <sep> is the separator in the output file.
        <replace_list> contains a double, list of items to be replaced,
            e.g. replace_list=[['.',',']] for European decimal notation.
        """
        self.ret_string = '' # string to be returned
        self.sep = sep
        self.replace_list = replace_list

        self.line_start = True    # a line is just starting

    def write_word(self, word):
        """
        Add a new entry to the csv.
        """
        if not self.line_start:
            self.ret_string += self.sep + word
        else:
            self.line_start = False
            self.ret_string += word

    def new_line(self):
        """
        Start a new line.
        """
        self.ret_string += '\n'
        self.line_start = True

    def write_line(self, words):
        """
        Writes a line with list <words>.
        """
        for word in words:
            self.write_word(word)
        self.new_line()

    def return_csv(self):
        """
        Return the table that has been written with write_line.
        """
        return self.ret_string

    def write_to_file(self, file_name):
        """
        Write the table to a file.
        """
        write_to_file(self.return_table(), file_name)
                      
class ky_file:
    """
    Handling of an input keystroke file.
    """
    def __init__(self,file_name):
        self.file_name = file_name
        try:
            self.fil = open(self.file_name, 'r')
        except IOError:
            self.fil = None
        self.out_list = []
        
    def def_input(self,text, default=''):
        if not self.fil == None:
            try:
                tmp = self.fil.next()[:-1]
            except StopIteration:
                tmp = default
        else:
            tmp = default
            
        ret_str = def_input(text,tmp)
        
        self.out_list.append(ret_str+'\n')
        return ret_str
        
    def close_write(self,wfile_name=None):
        """
        Overwrite the file with the new key strokes.
        """
        if wfile_name == None: wfile_name = self.file_name
            
        if not self.fil == None: self.fil.close()
        self.fil = open(wfile_name,'w')
        self.fil.writelines(self.out_list)
        self.fil.close()

def write_to_file(string, file_name):
    """
    Write <string> to file with <file_name>.
    """
    w_file = open(file_name, 'w')
    w_file.write(string)
    w_file.close()
    
def def_input(text, default=''):
    """
    Ask for user input.
    If specified a default is taken when the input is empty.
    """
    ask_str = text
    if default != '': ask_str += ' [%s]'%default
    ask_str += ':'
    
    print(ask_str)
    tmp = raw_input('')
    if tmp == '':
        tmp = default
        
    return tmp
        
if __name__ == '__main__':
    cm = csv_maker(sep=';')
    cm.write_line(['a', 'bea', 'da'])
    cm.write_word('3')
    cm.write_word('r')
    cm.new_line()
    cm.write_word('t')
    print(cm.return_csv())

    #cm.write_to_file('test.txt')

