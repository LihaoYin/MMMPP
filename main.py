import argparse

from Dataset import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_path", type=str, default="datasets.csv", help="Path to training dataset")
    parser.add_argument("--label_path", type=str, default="", help="Path to labels of training dataset")
    parser.add_argument("--test_dataset_path", type=str, default="test.csv", help="Path to testing dataset")
    parser.add_argument("--test_label_path", type=str, default="test.csv", help="Path to labels of testing dataset")
    parser.add_argument("--model", type=str, default="single_level", help="Type of model")
    parser.add_argument("--naacounts", type=int, default=1000, help="Number of accounts")
    parser.add_argument("--mdays", type=int, default=1, help="Number of days")
    parser.add_argument("--event_types", type=int, default=1, help="Number of event types")
    parser.add_argument("--time_slot", type=tuple, default=(0,1), help="Time interval for each repeats")

    opt = parser.parse_args()

    train_dataset = Dataset(
        dataset_path = opt.dataset_path,
        label_path = opt.label_path,
        model = opt.model,
        naccounts = opt.naccounts,
        mdays = opt.mdays
    )

    nclusters = 2
    train_dataset.ES_initial(nclusters)
    id, log_likelihood = train_dataset.MSMPP(nclusters)
