filename = r'C:\pythontest\Lib\site-packages\numpy\distutils\fcompiler\gnu.py'
with open(filename) as f:
     lines = f.read().replace('raise NotImplementedError("Only MS compiler supported with gfortran on win64")', 'pass')
with open(filename, "w") as f1:
     f1.write(lines)
