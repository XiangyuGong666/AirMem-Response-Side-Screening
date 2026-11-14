"""Desktop and command-line entry point."""
import argparse
from pathlib import Path


def main(argv=None):
    parser=argparse.ArgumentParser(description='AirMem Response-Side Screening · V19.0')
    sub=parser.add_subparsers(dest='command')
    scan=sub.add_parser('calculate', help='Analyze the 10 Hz acceleration and wind records')
    scan.add_argument('data', type=Path, help='Extracted 2022Q4 directory or its parent')
    scan.add_argument('output', type=Path, help='New or empty run directory')
    scan.add_argument('--original-time', type=Path, help='Optional timestamped MAT directory or time-vector ZIP for Figure 3')
    plot=sub.add_parser('plot', help='Generate selected English figures from a completed run')
    plot.add_argument('output', type=Path, help='Completed run directory')
    plot.add_argument('figure', nargs='+', help='Figure numbers, or all')
    plot.add_argument('--data', type=Path, help='Relocated 2022Q4 inputs, if necessary')
    args=parser.parse_args(argv)
    if args.command=='calculate':
        from .analysis import run
        print(run(args.data,args.output,original_root=args.original_time))
    elif args.command=='plot':
        from .figures import TITLES,build,export,read_run
        run=read_run(args.output)
        all_figures=args.figure==['all']
        if all_figures:
            numbers=list(TITLES)
            if 'coverage_daily.csv' not in run['outputs']:
                numbers.remove(3)
                print('Figure 3 skipped: no time index was supplied.')
        else:
            try: numbers=[int(n) for n in args.figure]
            except ValueError: plot.error('Supply figure numbers or the word all.')
        for n in numbers:
            fig=build(n,args.output,args.data)
            png,_=export(fig,n,args.output); print(png,flush=True)
    else:
        from .desktop import start
        return start()
    return 0


if __name__=='__main__':
    raise SystemExit(main())
