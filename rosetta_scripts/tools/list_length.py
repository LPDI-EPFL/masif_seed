import os,sys,string,gzip

def find_pdb_files():

    list_of_files=[]
    cmd = 'find ./results/*/* -name "*.pdb" -print '
    for file in os.popen(cmd).readlines():
        name = file[:-1]
        list_of_files.append(name)

    return list_of_files 


def extract_length( file ):

    list_resi=[]
    input_file = open(file)
    for line in input_file.readlines():
        split_line=line.split()

        if len(split_line) > 3 and split_line[0] == "ATOM" and split_line[4]=='B':

            list_resi.append(split_line[5])
    length = int(list_resi[-1]) - int(list_resi[0]) + 1
    return length


list_of_files = find_pdb_files()
i=0
with open('list_length.txt','w') as outfile:
    for filename in list_of_files:
        length=extract_length(filename)
        short=filename.split('/')[-1].strip('.pdb')
        output=short + '\t' + str(length) + '\n'
        outfile.write(output)
        i+=1
        print(str(i)+' out of '+str(len(list_of_files)))
