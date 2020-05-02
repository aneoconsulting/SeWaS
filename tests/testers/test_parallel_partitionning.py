import unittest

import subprocess
import re
import os

import adios2

class TestParallelPartitionning(unittest.TestCase):

    def setUp(self):
        try:
            self.mpiexec_executable = os.environ['MPIEXEC_EXECUTABLE']
            self.mpiexec_numproc_flag = os.environ['MPIEXEC_NUMPROC_FLAG']
        except:
            self.mpiexec_executable = "mpirun"
            self.mpiexec_numproc_flag = "-n"

    def test_parallel_partitionning_over_x_axis(self):

        for P in [2, 3, 4, 5]:
            test_result = subprocess.run(self.mpiexec_executable + " " + self.mpiexec_numproc_flag + " " + str(P) + " sewas --cx 40 --cy 200 --cz 100 --P " + str(P) + " --Q 1 --R 1 --nthreads 1 --dfile=tests/TestC.json",
                                         shell=True,
                                         capture_output=True,
                                         encoding='utf8')

            self.assertEqual(test_result.returncode, 0,
                             "[parallel partitionning over x axis] Failed to run test case")

            [vx] = [float(v) for v in re.findall(
                r'Vx\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]

            self.assertAlmostEqual(
                vx, 0.005981223684110059, msg="[parallel partitionning over x axis] failed for P = " + str(P))


if __name__ == '__main__':
    unittest.main()
