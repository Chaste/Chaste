import sys
from cxxtest import cxxtestgen 

def main(args=sys.argv):
    print("Running helper")
    cxxtestgen.main(args)

if __name__ == '__main__':
    main(sys.argv)