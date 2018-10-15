import os
import pytest
from csm.main.csm_run import csm_run
from tests.output_tests.conftest import test_folder
from tests.output_tests.utils import CheckFolder, standard_folder


class TestYaffa(CheckFolder):
    #@pytest.fixture(scope="class",
    #                params=[
    #                    "s6-fibonacci"
    #                ])
    #def test_name(self, request):
    #    return request.param

    def test_run(self, test_name):
        # run csm, saving results and checking if it crashes
        in_dir=os.path.join(test_folder, test_name)
        out_dir=self.result_folder(test_name)
        with open(os.path.join(in_dir, "testcommand.txt"), 'r') as file:
            command=file.read()
        os.chdir(in_dir)
        command=command.replace("__OUTPUT__", out_dir)
        result=csm_run(command.split())
        assert result
