import subprocess
import pkg_resources
import sys


def get_header(infile):
    """
    Read fragment file header
    """
    header = []
    with open(infile, 'r') as i:
        while True:
            l = i.readline()
            if l.startswith("#"):
                l = l.strip("\n")
                l = l.split("\t")
                kv = {}
                for tag in l:
                    if tag.startswith("#"):
                        kv['key'] = tag
                    else:
                        attributes = tag.split(":")
                        kv[attributes[0]] = attributes[1]
                header.append(kv)
            else:
                break
    return header


def update_header(header):
    """
    Add information to header
    """
    if len(header) == 0:
        return header # no header
    else:
        if header[0]['key'] == "#HD":
            header[0]['SO'] = 'coordinate'
    # add PG lines
    version = pkg_resources.require("sinto")[0].version
    pg = {
        'key': '#PG',
        'ID': 'sinto',
        'PN': 'sinto.sort',
        'VN': version,
        'CL': ' '.join(sys.argv)
    }
    # check if ID already taken
    counter = 1
    for i in header:
        if i['key'] == '#PG':
            if i['ID'] == 'sinto':
                i['ID'] = 'sinto.' + str(counter)
                i['PP'] = 'sinto'
                counter += 1
    header.append(pg)
    return header


def write_header(header, outfile):
    with open(outfile, 'w') if outfile != "-" else sys.stdout as o:
        for i in header:
            o.write(i['key'])
            for k in i.keys():
                if k == 'key':
                    pass
                else:
                    o.write("\t" + k + ":" + i[k])
            o.write("\n")


def sort_frags(infile, outfile):
    header = get_header(infile)
    header = update_header(header)
    write_header(header=header, outfile=outfile)
    grep_command = ["grep", "-v", '^#', infile]
    sort_command = ["sort", "-k1,1", "-k2,2n"]
    
    grep_ps = subprocess.Popen(grep_command, stdout=subprocess.PIPE)
    if outfile == "-":
        outp=sys.stdout
        subprocess.call(sort_command, stdin=grep_ps.stdout)
    else:
        outp = open(outfile, 'a')
        subprocess.call(sort_command, stdin=grep_ps.stdout, stdout=outp)
        outp.close()
