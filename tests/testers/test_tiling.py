import unittest

import subprocess
import re
import os

import sys

import adios2

class TestTiling(unittest.TestCase):
    
    def setUp(self):
        self.test_env = os.environ['PYTHONPATH']
        print(self.test_env)

    def test_tiling_over_x_axis(self):
        for cx in range(3, 100, 5):
            test_result = subprocess.run("sewas --cx " + str(cx) + " --cy 200 --cz 100 --P 1 --Q 1 --R 1 --nthreads 2 --dfile=tests/TestC.json",
                                         shell=True,
                                         capture_output=True,
                                         encoding='utf8')

            self.assertEqual(test_result.returncode, 0,
                             "[tiling over x axis] Failed to run test case")

            [vx] = [float(v) for v in re.findall(
                r'Vx\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vy] = [float(v) for v in re.findall(
                r'Vy\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vz] = [float(v) for v in re.findall(
                r'Vz\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]

            self.assertAlmostEqual(
                vx, 0.005981223684110059, msg="[tiling over x axis] failed for cx = " + str(cx))
            self.assertAlmostEqual(
                vy, 0.005981223684110059, msg="[tiling over x axis] failed for cx = " + str(cx))
            self.assertAlmostEqual(
                vz, 0.0009125988001892094, msg="[tiling over x axis] failed for cx = " + str(cx))

    def test_tiling_over_y_axis(self):
        for cy in range(3, 100, 5):
            test_result = subprocess.run("sewas --cx 200 --cy " + str(cy) + " --cz 100 --P 1 --Q 1 --R 1 --nthreads 2 --dfile=tests/TestC.json",
                                         shell=True,
                                         capture_output=True,
                                         encoding='utf8')

            self.assertEqual(test_result.returncode, 0,
                             "[tiling over y axis] Failed to run test case")

            [vx] = [float(v) for v in re.findall(
                r'Vx\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vy] = [float(v) for v in re.findall(
                r'Vy\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vz] = [float(v) for v in re.findall(
                r'Vz\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]

            self.assertAlmostEqual(
                vx, 0.005981223684110059, msg="[tiling over y axis] failed for cy = " + str(cy))
            self.assertAlmostEqual(
                vy, 0.005981223684110059, msg="[tiling over y axis] failed for cy = " + str(cy))
            self.assertAlmostEqual(
                vz, 0.0009125988001892094, msg="[tiling over y axis] failed for cy = " + str(cy))

    def test_tiling_over_z_axis(self):
        for cz in range(3, 50, 5):
            test_result = subprocess.run("sewas --cx 200 --cy 200 --cz " + str(cz) + " --P 1 --Q 1 --R 1 --nthreads 2 --dfile=tests/TestC.json",
                                         shell=True,
                                         capture_output=True,
                                         encoding='utf8')

            self.assertEqual(test_result.returncode, 0,
                             "[tiling over z axis] Failed to run test case")

            [vx] = [float(v) for v in re.findall(
                r'Vx\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vy] = [float(v) for v in re.findall(
                r'Vy\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]
            [vz] = [float(v) for v in re.findall(
                r'Vz\|\| = ([0-9]\.[0-9]+[Ee]?[+-]?[0-9]*)', test_result.stdout)]

            self.assertAlmostEqual(
                vx, 0.005981223684110059, msg="[tiling over z axis] failed for cx = " + str(cz))
            self.assertAlmostEqual(
                vy, 0.005981223684110059, msg="[tiling over z axis] failed for cx = " + str(cz))
            self.assertAlmostEqual(
                vz, 0.0009125988001892094, msg="[tiling over z axis] failed for cx = " + str(cz))


if __name__ == '__main__':
    unittest.main()
