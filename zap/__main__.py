# -*- coding: utf-8 -*-

import argparse
import logging
import sys

from .version import __version__
from .zap import process


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='ZAP (the Zurich Atmosphere Purge) is a high precision sky'
                    ' subtraction tool.'
    )
    parser.add_argument('incube', help='Input datacube path')
    parser.add_argument('--version', '-V', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('--debug', '-d', action='store_true',
                        help='Show debug info')
    parser.add_argument('--no-clean', action='store_true',
                        help='Disable NaN values interpolation')
    parser.add_argument('--outcube', '-o', default='DATACUBE_FINAL_ZAP.fits',
                        help='Output datacube path')
    parser.add_argument('--skycube', help='Sky datacube path')
    parser.add_argument('--cfwidth', type=int, default=100,
                        help='Window size for the continuum filter')
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        process(args.incube, outcubefits=args.outcube, clean=not args.no_clean,
                skycubefits=args.skycube, cfwidth=args.cfwidth)
    except KeyboardInterrupt:
        sys.exit('Interrupted!')
    except Exception as e:
        sys.exit('Failed to process file: %s' % e)


if __name__ == "__main__":
    main()
