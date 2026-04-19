"""
Pipeline execution and operation queuing.

Provides Pipeline class for chaining data processing operations into
reproducible workflows with automatic logging and progress tracking.
"""

import copy
import datetime
import gc
import os
import time

import numpy as np
import pandas as pd

from clibas.baseclasses import Handler


class Pipeline(Handler):
    """
    Operation queue and execution manager for data processing workflows.

    Provides infrastructure for building, executing, and monitoring data
    processing pipelines. Operations are enqueued and executed sequentially
    with automatic timing, logging, and progress tracking.

    Example:
        >>> import clibas as C
        >>> C.initialize_from_config('config.yaml')
        >>>
        >>> #build pipeline
        >>> C.pipeline.enque([
        ...                   C.fastq_parser.translate(),
        ...                   C.fastq_parser.len_filter(where='pep'),
        ...                   C.fastq_parser.save(where='pep', fmt='fasta')
        ... ])
        >>>
        >>> #execute
        >>> loader = C.data_loader.fetch_gz_from_dir(data_dir='./sequencing_data/')
        >>> data = C.pipeline.load_and_run(loader=loader)
    """

    def __init__(self, *args):
        super(Pipeline, self).__init__(*args)
        self._on_startup()
        return

    def __repr__(self):
        return f"<Pipeline object; current queue size: {len(self.que)} op(s)>"

    def _on_startup(self):
        self.que = []
        if not hasattr(self, "exp_name"):
            self.exp_name = "untitled_exp"
        return

    def _describe_data(self, data=None):
        """
        Go over every dataset for every sample and log all array shapes.
        Used during dequeing to keep track of the data flow."""
        data_descr = []

        if data is None:
            return data_descr

        for sample in data:
            data_descr.append((sample.name, sample.size))
            if sample.size:
                msg = f"{sample.name} dataset size: {sample.size}"
                self.logger.info(msg)
            else:
                msg = 65 * "-"
                msg = f"Sample {sample.name} has zero entries remaining!"
                self.logger.warning(msg)

        msg = 65 * "-"
        self.logger.info(msg)

        return data_descr

    def _reassemble_summary(self, summary):
        ops = []
        times = []
        samples = []

        # the code below is a mess, but the task is trivial,
        # so whatever; fix if nothing better to do
        for x in summary:
            ops.append(x["op"])
            times.append(x["op_time"])
            for j in x["data_description"]:
                samples.append(j[0])

        samples = list(set(samples))
        sizes = np.zeros((len(summary), len(samples)))
        for i, entry in enumerate(summary):
            for tup in entry["data_description"]:
                for j, name in enumerate(samples):
                    if tup[0] == name:
                        sizes[i, j] = tup[1]

        df = pd.DataFrame(columns=["elapsed time, s"] + samples, index=ops)
        df["elapsed time, s"] = times
        for i, name in enumerate(samples):
            df[name] = sizes[:, i]

        return df

    def enque(self, ops):
        """
        Add operations to the pipeline queue.

        Args:
            ops (list): List of callable operations. Each operation should
                accept a Data object and return a transformed Data object.

        Example:
            >>> C.pipeline.enque([
            ...                   parser.translate(force_at_frame=0),
            ...                   parser.len_filter(where='pep'),
            ...                   parser.save(where='pep', fmt='fasta')
            ... ])
        """
        for func in ops:
            self.que.append(func)

        msg = (
            f"{len(ops)} ops appended to pipeline; current queue size: {len(self.que)}"
        )
        self.logger.info(msg)
        return

    def run(self, data=None, save_summary=True):
        """
        Execute the enqueued pipeline.

        Processes operations sequentially, transforming data through each step.
        Logs timing and data size information at each stage.

        Args:
            data (Data, optional): Input data. If None, first operation must
                load data (e.g., a loader function).

            save_summary (bool): If True, saves CSV summary of pipeline
                execution to logs directory. Default is True.

        Returns:
            Data: Transformed Data object after finishing all operations.

        Example:
            >>> data = C.pipeline.run()  #execute full pipeline
        """
        if len(self.que) == 0:
            msg = "<Pipeline>: the queue is empty! Nothing to run."
            self.logger.warning(msg)
            return

        summary = list()
        data_descr = self._describe_data(data)
        summary.append({"op": None, "op_time": None, "data_description": data_descr})

        for _ in range(len(self.que)):
            func = self.que.pop(0)
            msg = f"Queuing <{func.__name__}> op. . ."
            self.logger.info(msg)

            t = time.time()
            data = func(data)
            op_time = np.round(time.time() - t, decimals=3)

            msg = f"The operation took {op_time} s"
            self.logger.info(msg)
            data_descr = self._describe_data(data)

            summary.append(
                {
                    "op": func.__name__,
                    "op_time": op_time,
                    "data_description": data_descr,
                }
            )

        if save_summary:
            summary = self._reassemble_summary(summary)
            n = datetime.datetime.now()
            timestamp = f"_{n.year}_{n.month}_{n.day}_{n.hour}_{n.minute}_{n.second}"

            fname = f"{self.exp_name}_pipeline_summary_{timestamp}.csv"
            self._prepare_destinations(root=self.dirs.logs)
            summary.to_csv(os.path.join(self.dirs.logs, fname))

        return data

    def load_and_run(self, loader, save_summary=True):
        """
        Execute pipeline with data loading as first step.

        Convenience method that prepends a loader function to the queue
        and executes the full pipeline.

        Args:
            loader (callable): Function that returns a Data object.

            save_summary (bool): If True, saves execution summary.

        Returns:
            Data: Transformed Data object after finishing all operations.

        Example:
            >>> data = C.pipeline.load_and_run(C.data_loader.fetch_gz_from_dir())
        """
        self.que = [
            loader,
        ] + self.que
        data = self.run(data=None, save_summary=save_summary)
        return data

    def stream(self, streamer, save_summary=True):
        """
        Execute pipeline on streamed data for memory efficiency.

        Processes samples one at a time from a generator, useful when total
        dataset exceeds available memory. Each sample processed independently
        through the full pipeline.

        Args:
            streamer (generator): Generator yielding Sample objects.

            save_summary (bool): If True, saves execution summary for each
                sample. Default is True.

        Example:
            >>> #process large dataset sample-by-sample
            >>> C.pipeline.stream(C.data_loader.stream_from_gz_dir())
        """
        from clibas.datatypes import Data

        que = copy.deepcopy(self.que)
        for sample in streamer:
            # setting the exp name will help writing
            # pipeline summaries for each sample
            self.exp_name = sample.name

            # turn sample into a Data instance and pass it though the pipeline
            data = Data([sample])
            self.que = copy.deepcopy(que)
            self.run(data=data, save_summary=save_summary)

            # unless the data is deleted, two datasets are stored in memory
            # before the next call to data when a new sample is made.
            del data
            gc.collect()

            # reset library designs back to the original
            if hasattr(self, "P_design"):
                self.P_design.rebuild()

            if hasattr(self, "D_design"):
                self.D_design.rebuild()

        # unset the experiment name after all done
        self.exp_name = None
        return
