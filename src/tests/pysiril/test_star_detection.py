from pysiril.siril import Siril
from pysiril.wrapper import Wrapper
import os, sys

def test_star_detection(path_to_siril = None):
    """ test number of stars and use of detection parameters """
    if path_to_siril is not None:
        app = Siril(siril_exe = path_to_siril)
    else:
        app = Siril()

    app.Open()
    cmd=Wrapper(app)

    cmd.set16bits()
    cmd.setext('fit')
    cmd.cd(os.getcwd())
    cmd.load('star_test1')
    cmd.findstar(out='star_list1.lst')

    # start of the star list file should be like that:
    line1 = "# 52 stars found using the following parameters:\n"
    line2 = "# sigma=1.00 roundness=0.50 radius=5 relax=0 profile=0 minbeta=1.5 max_stars=0 layer=0 minA=0.00 maxA=0.00 maxR=1.00\n"

    with open('star_list1.lst') as file:
        line = file.readline()
        assert line1 == line, f"line is {line}"
        line = file.readline()
        assert line2 == line, f"line is {line}"

    # retrying with different parameters
    sigma = 0.7
    roundness = 0.75
    convergence = 3
    app.Execute(f"setfindstar -sigma={sigma} -roundness={roundness} -convergence={convergence}")
    cmd.findstar(out='star_list2.lst')

    line1 = "# 62 stars found using the following parameters:\n"
    line2 = "# sigma=0.70 roundness=0.75 radius=5 relax=0 profile=0 minbeta=1.5 max_stars=0 layer=0 minA=0.00 maxA=0.00 maxR=1.00\n"
    with open('star_list2.lst') as file:
        line = file.readline()
        assert line1 == line, f"line is {line}"
        line = file.readline()
        assert line2 == line, f"line is {line}"

    app.Close()
    del app

if __name__ == '__main__':
    os.chdir(os.path.dirname(sys.argv[0]))
    if len(sys.argv) > 1:
        test_star_detection(sys.argv[1])
    else:
        test_star_detection()

