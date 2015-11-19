import fileinput
file = r'C:\pythontest\Lib\site-packages\numpy\distutils\fcompiler\gnu.py'
for line in fileinput.FileInput(file, inplace=True):
    search = 'raise NotImplementedError("Only MS compiler supported with gfortran on win64")'
    if search in line:
        line = line.replace(search, 'pass')
        print('replaced', search, 'in', file)
