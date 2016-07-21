#!/usr/bin/env python2

# Script created by Bobby Kim and Nick Schafer on 1/8/13
# This script assumes that you have a file called "pdbtmall"
# in the current directory. This file is an XML version of the
# PDBTM database, which can be found here:
# http://pdbtm.enzim.hu/?_=/download/files
# The script takes two arguments:
# python getTransmembraneRegions.py pdbidlist
# It then searches for transmembrane regions that correspond to
# that pdbids in the pdbidlist and writes the information out to a
# file called "pdbidCHAINID.tm". The format of that file is a
# single column of 1s and 0s, each row corresponding to one
# residue, where 1 means that it is transmembrane and 0 means
# that it is not. WARNING: In its original form, the script
# searches for transmembrane HELICES only, but this can be
# changed by searching for other "regiontype"s. For information
# about other region types, see:
# http://pdbtm.enzim.hu/?_=/docs/manual

import sys
import xml.etree.ElementTree as ET
import os.path

def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]

def cyto_peri_switcher():
    global cyto_peri_switch
    if cyto_peri_switch == 1:
        cyto_peri_switch = 3
    elif cyto_peri_switch == 3:
        cyto_peri_switch = 1

cyto_peri_switch = 3
# import data from file
pdbfile = open("pdbtmall","r")
pdbtm = ET.parse(pdbfile)
remove_namespace(pdbtm, 'http://pdbtm.enzim.hu')

idlistfile=sys.argv[1]
idlist=open(idlistfile, "r")
for line in iter(idlist):
    # find transmembrane regions and print them out
    print line
    pdbid=line[0:4]
    print pdbid
    chainid=line[5]
    print chainid
    outputfilename='./%s%s.tm' % (pdbid, chainid)

    if os.path.isfile(outputfilename):
        print "%s already exists, skipping" % outputfilename
        continue

    outputfile=open(outputfilename,"w")
    prev_regiontype = ''
    transmembraneregions=pdbtm.findall("./pdbtm[@ID='%s']/CHAIN[@CHAINID='%s']/REGION" % (pdbid, chainid))
    print pdbid, chainid
    for region in transmembraneregions:
        start = int(region.get('seq_beg'))
        end = int(region.get('seq_end'))
        regiontype = region.get('type')
        print region.get('seq_beg'), region.get('seq_end'), region.get('type')
        for i in range(start,end+1):
            if regiontype == 'H':
                outputfile.write('2\n')
            else:
                if prev_regiontype == 'H':
                    cyto_peri_switcher()
                if cyto_peri_switch == 1:
                    outputfile.write('1\n')
                elif cyto_peri_switch == 3:
                    outputfile.write('3\n')
            prev_regiontype = regiontype
    outputfile.close()
