import argparse

from clipseq_pipeline.analysis import run_analysis
from clipseq_pipeline.prediction import run_prediction
from clipseq_pipeline.processing import run_processing

STAGES = {
    "processing": run_processing,
    "prediction": run_prediction,
    "analysis": run_analysis
}

def main():
    parser = argparse.ArgumentParser(description="Run pipeline stages.")
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=STAGES.keys(),
        default=list(STAGES.keys()),
        help="Stages to run (default: all)"
    )
    args = parser.parse_args()

    for stage in args.stages:
        print(f"\n=== Running stage: {stage} ===")
        STAGES[stage]()
        print(f"=== Finished stage: {stage} ===")

if __name__ == "__main__":
    main()