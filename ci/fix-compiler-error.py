# Handles GCC and CYGWIN Compilers
import os

filename = {
     'gcc': r'C:\pythontest\Lib\site-packages\numpy\distutils\fcompiler\gnu.py',
     'cygwin': r'C:\pythontest\Lib\distutils\cygwinccompiler.py'}

err_msg = {'gcc': ''.join(['raise NotImplementedError("Only MS compiler ',
                           'supported with gfortran on win64")']),
           'cygwin': "            return ['msvcr100']"}
new_line = {'gcc': 'pass',
            'cygwin': '\n'.join(["            return ['msvcr100']",
                                 "        elif msc_ver == '1900':",
                                 "            return ['vcruntime140']"])}
for ckey in filename.keys():
    if os.path.isfile(filename[ckey]):
        # Open the desired file and replace bad lines, if they exist
        with open(filename[ckey], 'r') as fin:
            lines = fin.read().replace(err_msg[ckey], new_line[ckey])

        # Overwrite file with modified text
        with open(filename[ckey], "w") as fout:
            fout.write(lines)
