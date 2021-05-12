# from context import ridkit
import glob, os
# ridkit.gen_rid.gen_rid('./tmp_rid', 'benchmark_mol', 'benchmark_json/rid_not_enough.json')
# ridkit.make_enhc.make_enhc(iter_index=0, 
#                json_file='/home/dongdong/wyz/ridkit-master/tests/benchmark_json/rid.json',
#                graph_files=[],
#                mol_dir='/home/dongdong/wyz/ridkit-master/tests/benchmark_mol',
#                cv_file='/home/dongdong/wyz/ridkit-master/tests/benchmark_json/cv.json',
#                base_dir='/home/dongdong/wyz/ridkit-master/tests/test_case/case_make_enhc')


# cwd = os.getcwd()
# os.chdir('/home/dongdong/wyz/refinement3_6w_320/R0968s2.run06/iter.000000/02.train/')
# file_list = glob.glob('*.pb')
os.chdir("/home/dongdong/wyz/ridkit-master/tests/test_case/case_make_enhc_iter1/iter.000000/02.train/")

for ff in file_list:
    os.system("touch {}".format(ff))