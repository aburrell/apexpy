filename = r'C:\pythontest\Lib\site-packages\numpy\distutils\fcompiler\gnu.py'
with open(filename) as f:
     lines = f.read().replace('raise NotImplementedError("Only MS compiler supported with gfortran on win64")', 'pass')
with open(filename, "w") as f1:
     f1.write(lines)

filename = r'C:\pythontest\Lib\distutils\cygwinccompiler.py'     
with open(filename) as f:
     lines = f.read().replace("            return ['msvcr100']", "            return ['msvcr100']\n        elif msc_ver == '1600':\n            return ['vcruntime140']")
     print(lines)
with open(filename, "w") as f1:
     f1.write(lines)
