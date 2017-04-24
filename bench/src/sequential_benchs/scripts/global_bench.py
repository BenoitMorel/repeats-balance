import sys, string, os 
import subprocess

current_lib_path = '../lib/current'
tree404 = '../../data/404/unrooted.newick'
seq404 = '../../data/404/404.phy'


class Result:
    full_traversal_dataset = []
    full_traversal_times = []

class InputLib:
    lib=''
    use_repeats='0'
    update_repeats='0'
    result = Result()

    def __init__(self, lib, use_repeats, update_repeats):
        self.lib = lib
        self.use_repeats = use_repeats
        self.update_repeats = update_repeats




def bench_lib(lib) :
    os.system('cp ' + lib.lib + '/* ' + current_lib_path)
    #os.system('make clean && make')
    
    time = subprocess.check_output(['./main', 'full_traversal',
        tree404, seq404, lib.use_repeats, lib.update_repeats, 
        '50', 'cpu'])
    print time




def output_latex(latex_file, libs):
    print 'yoyo'


# BEGIN

lib1 = InputLib('../lib/libpll_sr_standard', '1', '0')
lib2 = InputLib('../lib/libpll_sr_standard', '0', '0')
input_libs = [lib1, lib2]



os.environ['LD_LIBRARY_PATH'] = current_lib_path
for lib in input_libs:
    bench_lib(lib)
output_latex('plop.tex', input_libs) 

print 'done'

