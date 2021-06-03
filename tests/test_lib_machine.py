from context import rid
import unittest, os, shutil, json, glob, filecmp
from rid.lib.machine import set_machine, set_resource
import mock

class TestUtilsMachine(unittest.TestCase):

    def setUp(self):
        self.case_dir = "test_case/case_utils"
        self.machine_json = os.path.join(self.case_dir, 'machine.json')

    @mock.patch('rid.lib.machine.Resources.load_from_dict')
    def test_set_resource(self, mock_resources):
        args = {
            "number_node": 1,
            "cpu_per_node": 8,
            "gpu_per_node": 1,
            "queue_name": 'V100',
            "group_size": 2,
            "if_cuda_multi_devices": False
        }
        set_resource(self.machine_json, target="abcd")
        mock_resources.assert_called_with(args)

    def test_set_machine(self):
        machine = set_machine(self.machine_json, target="pbs")
        self.assertTrue(type(machine).__name__ == "PBS")
        machine = set_machine(self.machine_json, target="slurm")
        self.assertTrue(type(machine).__name__ == "Slurm")


if __name__ == "__main__":
    unittest.main()