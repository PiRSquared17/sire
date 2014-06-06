
from Sire.Base import *


def test_ls(verbose = False):
    if verbose:
        print("Running process \"ls\"")

    p = Process.run("ls")
    p.wait()

    if verbose:
        print("command was successful!")

    if verbose:
        print("Running process \"echo 'process ran successfully'\"")

    p = Process.run("echo", "process ran successfully")
    p.wait()

if __name__ == "__main__":
    test_ls(True)

