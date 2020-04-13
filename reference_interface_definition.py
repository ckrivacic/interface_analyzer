'''
Define reference interface. 
Required: reference_interfaces  (list of list), 
reference_pdb (path to reference pdb file)
'''
import os

reference_interface_all = [126,127,125,19,124,80,79,128,78,18,76,129,11,16,
        14, 81, 13, 116, 12, 115, 15, 9, 8, 7]

reference_interfaces = [[],[]]
for resnum in reference_interface_all:
    if resnum < 90:
        reference_interfaces[0].append(resnum)
    else:
        reference_interfaces[1].append(resnum)

reference_interfaces.append(reference_interface_all)

reference_pdb = os.path.join('test_inputs', 'reference.pdb')


if __name__=="__main__":
    from utils import *
    for i in reference_interfaces:
        print(reslist_to_selstr(i))
