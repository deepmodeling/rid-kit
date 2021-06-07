from context import rid
import unittest, os, shutil, json, glob, filecmp
from rid.lib import utils
import mock


file_content = """this is a test
dpfe: DEEPFE TRUST_LVL_1=2.000000 TRUST_LVL_2=3.000000 MODEL= ARG=dih-000-01
MODEL_TEST=456 """
file_content_replaced = """this is a test
dpfe: DEEPFE TRUST_LVL_1=2.000000 TRUST_LVL_2=100 MODEL=test.pd ARG=dih-000-01
MODEL_TEST=456 """

class TestUtilsString(unittest.TestCase):
    def setUp(self):
        self.case_dir = "test_case/case_utils"
        pass
        
    def test_make_name(self):
        self.assertTrue((utils.make_iter_name(5) == "iter.000005"))
        self.assertTrue((utils.make_walker_name(5) == "005"))
        pass
    
    def test_repeat_to_length(self):
        astring = utils.repeat_to_length("rid_is_powerful_", 3)
        self.assertTrue(astring == "rid_is_powerful_rid_is_powerful_rid_is_powerful_")

    def test_cmd_append_log(self):
        cmd_log = utils.cmd_append_log("gmx mdrun -deffnm em", 'md.log')
        _cmd_log = "gmx mdrun -deffnm em 1> md.log 2> md.log"
        self.assertTrue(cmd_log, _cmd_log)
    
    def test_print_list(self):
        astring = utils.print_list(['test_rid', 'test_utils', 'test_print'], suffix='.txt')
        self.assertTrue(astring, "test_rid.txt,test_utils.txt,test_print.txt")

    def test_print_repeat(self):
        astring = utils.print_repeat_list(3, "test_utils")
        self.assertTrue(astring, "test_utils,test_utils,test_utils")


class TestUtilsFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.path.exists("tmp_utils/"):
            shutil.rmtree("tmp_utils")
        os.mkdir("tmp_utils")
        if not os.path.exists("test_case/case_utils/record.txt"):
            fp = open("test_case/case_utils/record.txt", "w")
            fp.write("0 0\n7 3\n")
            fp.close()
        if not os.path.exists("test_case/case_utils/rid_foo.txt"):
            fp = open("test_case/case_utils/rid_foo.txt", "w")
            fp.write("1")
            fp.close()
        if not os.path.exists("test_case/case_utils/rid_foo2.txt"):
            fp = open("test_case/case_utils/rid_foo2.txt", "w")
            fp.write("")
            fp.close()


    def setUp(self):
        self.test_dir = "tmp_utils/"
        self.case_dir = "test_case/case_utils"
        self.cwd = os.getcwd()
        pass

    def test_replace(self):
        with open(self.test_dir+"replace.txt", 'w') as test_file:
            test_file.write(file_content)
        utils.replace(self.test_dir+"replace.txt", 'MODEL=[^ ]* ', "MODEL=test.pd ")
        utils.replace(self.test_dir+"replace.txt", 'TRUST_LVL_2=[^ ]* ', "TRUST_LVL_2=100 ")
        with open(self.test_dir+"replace.txt", 'r') as test_file:
            ret = test_file.read()
        self.assertTrue(ret == file_content_replaced )
        pass

    def test_creat_path(self):
        _test_path = os.path.join(self.test_dir + "test_foo")
        utils.create_path(_test_path)
        self.assertTrue(os.path.exists(_test_path))

        _test_path_existed = os.path.join(self.test_dir + "test_existed")
        if os.path.exists(_test_path_existed):
            shutil.rmtree(_test_path_existed)
        os.mkdir(_test_path_existed)
        fp = open(_test_path_existed + "/test_rid.txt", 'w').close()
        if not os.path.exists(_test_path_existed + ".bk000"):
            os.mkdir(_test_path_existed + ".bk000")
        utils.create_path(_test_path_existed)
        self.assertTrue(os.path.exists(_test_path_existed))
        self.assertTrue(os.path.exists(_test_path_existed + ".bk001/test_rid.txt"))
        pass

    def test_copy_file_list(self):
        from_path = 'test_case/case_utils'
        file_list = ['rid_foo.txt', 'sub_dir']
        to_path = self.test_dir
        utils.copy_file_list(file_list, from_path, to_path)
        os.chdir(to_path)
        self.assertTrue(os.path.isfile('rid_foo.txt'))
        self.assertTrue(os.path.isfile('sub_dir/test_file_2.txt'))
        self.assertTrue(os.path.isfile('sub_dir/test_file.txt'))
        os.chdir(self.cwd)

    def test_checkfile(self):
        print("\n ! NOTE that 'checkfile' function has not been checked.\n")
        pass

    def test_record_iter(self):
        record_file = os.path.join(self.case_dir, "record.txt")
        with open(record_file, 'w') as record:
            record.write("0 0\n")
        utils.record_iter(record_file, 7, 3)
        with open(record_file, 'r') as record:
            ret = record.read()
        self.assertTrue(ret == "0 0\n7 3\n")
        pass

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("tmp_utils/"):
            shutil.rmtree("tmp_utils")
        if os.path.exists("test_case/case_utils/record.txt"):
            os.remove("test_case/case_utils/record.txt")
        if os.path.exists("test_case/case_utils/rid_foo.txt"):
            os.remove("test_case/case_utils/rid_foo.txt")
        if os.path.exists("test_case/case_utils/rid_foo2.txt"):
            os.remove("test_case/case_utils/rid_foo2.txt")




if __name__ == "__main__":
    unittest.main()