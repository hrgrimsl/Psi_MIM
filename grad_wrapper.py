import quantum_methods
import frag_script
import argparse, re
import xml.etree.ElementTree as ET


parser = argparse.ArgumentParser(description='Creates a set of fragments from a cml file',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('name', nargs=1, help='CML file to process')
parser.add_argument('-nm','--etam', type=int, default=1, help='Small eta for retrieving correlation energy?  (m for mini)', required=False)
parser.add_argument('-nw','--etaw', type=int, default=100, help='Large eta for retrieving hf energy? (w for wumbo)', required=False)
parser.add_argument('-tm', '--theorym', type = str, default = 'ccsd', help = 'High level theory method', required = False)
parser.add_argument('-tw', '--theoryw', type = str, default = 'scf', help = 'Low level theory method', required = False)
parser.add_argument('-b', '--basis', type = str, default = '6-31G*', help = 'basis', required = False)
parser.add_argument('-s','--scratch', type=str, help='Name of scratch directory', default="./scr", required=False)
args = vars(parser.parse_args())
small_good_args = {
'basis' : args['basis'],
'name' : args['name'],
'eta' : args['etam'],
'method' : args['theorym'],
'scratch': args['scratch']+"/smallgood"
}
small_bad_args = {
'basis' : args['basis'],
'name' : args['name'],
'eta' : args['etam'],
'method' : args['theoryw'],
'scratch': args['scratch']+"/smallbad"
}
big_args = {
'basis' : args['basis'],
'name' : args['name'],
'eta' : args['etaw'],
'method' : args['theoryw'],
'scratch': args['scratch']+"/big"
}
file_name = args['name'][0]
tree = ET.parse(file_name)
root = tree.getroot()

frag_script.Fragment(big_args)
#energy = quantum_methods.Compute_Energy(big_args)
#gradient = quantum_methods.Compute_Gradient(big_args, root)
frag_script.Fragment(small_bad_args)
#energy2 = -quantum_methods.Compute_Energy(small_bad_args)
#gradient2 = -1*quantum_methods.Compute_Gradient(small_bad_args, root)
frag_script.Fragment(small_good_args)
#energy3 = quantum_methods.Compute_Energy(small_good_args)
#gradient3 = quantum_methods.Compute_Gradient(small_good_args, root)
eg = quantum_methods.Compute_Gradient(args, root)
true_energy = eg[1]
true_gradient = eg[0]
#true_energy = energy+energy2+energy3
print("True Grad:")
print(true_gradient)
write_file = open(args['scratch']+"/"+'timshel', 'w')
#write_file.write(energy)
for i in range(0, len(true_gradient)):
    write_file.write(str(true_gradient[i][0])+" "+str(true_gradient[i][1])+" "+str(true_gradient[i][2])+"\n")

write_file_2 = open(args['scratch']+"/"+'amen', 'w')
print(true_energy)
write_file_2.write(str(true_energy))
