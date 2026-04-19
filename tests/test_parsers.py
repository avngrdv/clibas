"""
Created on Sat Oct 18 01:38:54 2025
@author: Alex Vinogradov
"""
from pathlib import Path


def test_pipeline_Heinis(Heinis_config, tmp_path, data_dir):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    C.pipeline.enque(
        [
            C.fastq_parser.trim_reads(left="CCAGCC", right="GGTTCTGGC", tol=1),
            C.fastq_parser.translate(stop_readthrough=False),
            C.analysis_tools.length_analysis(where="pep"),
            C.analysis_tools.length_analysis(where="dna"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.len_filter(where="pep"),
            C.fastq_parser.cr_filter(where="pep", loc=[0, 2], tol=2),
            C.fastq_parser.vr_filter(where="pep", loc=[1], sets=[1, 2]),
            C.analysis_tools.q_score_analysis(loc=None),
            C.fastq_parser.q_score_filt(minQ=30, loc=[1]),
            C.fastq_parser.fetch_at(where="pep", loc=[1]),
            C.fastq_parser.filt_ambiguous(where="pep"),
            C.analysis_tools.sequence_convergence_analysis(where="pep"),
            C.analysis_tools.token_convergence_analysis(where="pep", loc=None),
            C.analysis_tools.token_convergence_analysis(where="dna", loc=None),
            C.fastq_parser.unpad(),
            C.fastq_parser.save(where="pep", fmt="npy"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="fasta"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="fasta"),
            C.fastq_parser.library_design_match(where="pep"),
            C.fastq_parser.dataset_wide_count_summary(where="pep", top_n=2000),
        ]
    )

    loader = C.data_loader.fetch_fastq_from_dir(data_dir=data_dir)
    data = C.pipeline.load_and_run(loader=loader, save_summary=True)

    assert data.size == 2
    assert data[0].size == 1414
    assert data[1].size == 1653
    return


def test_pipeline_Walport(Walport_config, tmp_path, data_dir):
    import clibas as C

    C.initialize(Walport_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    C.pipeline.enque(
        [
            C.fastq_parser.translate(stop_readthrough=False),
            C.analysis_tools.length_analysis(where="pep"),
            C.analysis_tools.length_analysis(where="dna"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.len_filter(where="pep"),
            C.fastq_parser.cr_filter(where="pep", loc=[1], tol=3),
            C.fastq_parser.vr_filter(where="pep", loc=[0], sets=[1, 2, 3]),
            C.analysis_tools.q_score_analysis(loc=[0]),
            C.fastq_parser.fetch_at(where="pep", loc=[0]),
            C.fastq_parser.filt_ambiguous(where="pep"),
            C.analysis_tools.sequence_convergence_analysis(where="pep"),
            C.analysis_tools.token_convergence_analysis(where="pep", loc=None),
            C.analysis_tools.token_convergence_analysis(where="dna", loc=None),
            C.fastq_parser.unpad(),
            C.fastq_parser.save(where="pep", fmt="npy"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="fasta"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="fasta"),
            C.fastq_parser.library_design_match(where="pep"),
            C.fastq_parser.dataset_wide_count_summary(where="pep", top_n=2000),
            C.fastq_parser.drop_data(where="Q"),
        ]
    )

    loader = C.data_loader.fetch_gz_from_dir(data_dir=data_dir)
    data = C.pipeline.load_and_run(loader=loader, save_summary=True)

    assert data.size == 2
    assert data[0].size == 1420
    assert data[1].size == 513

def test_no_config_parser(tmp_path, data_dir):
    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "outputs")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "outputs")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "outputs")

    (tmp_path / "outputs").mkdir()

    C.pipeline.enque(
        [
            C.fastq_parser.trim_reads(left="CCAGCC", right="GGTTCTGGC", tol=1),
            C.fastq_parser.translate(),
            C.analysis_tools.length_analysis(where="pep"),
            C.analysis_tools.length_analysis(where="dna"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.analysis_tools.q_score_analysis(loc=None),
            C.fastq_parser.q_score_filt(avgQ=55, loc=None),
            C.fastq_parser.q_score_filt(minQ=10, loc=None),
            C.fastq_parser.filt_ambiguous(where="pep"),
            C.analysis_tools.sequence_convergence_analysis(where="pep"),
            C.analysis_tools.token_convergence_analysis(where="pep", loc=None),
            C.analysis_tools.token_convergence_analysis(where="dna", loc=None),
            C.fastq_parser.unpad(),
            C.fastq_parser.save(where="pep", fmt="npy"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="fasta"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="csv"),
            C.fastq_parser.count_summary(where="dna", top_n=1000, fmt="fasta"),
            C.fastq_parser.dataset_wide_count_summary(where="pep", top_n=2000),
        ]
    )
    
    loader = C.data_loader.fetch_fastq_from_dir(data_dir=data_dir)
    data = C.pipeline.load_and_run(loader=loader, save_summary=True)
    
    assert data.size == 2
    assert data[0].size == 136
    assert data[1].size == 1859
    
    
def test_dir_streaming(Heinis_config, tmp_path, data_dir):
    import clibas as C

    C.initialize(Heinis_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "streaming_dir")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "streaming_dir")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "streaming_dir")

    (tmp_path / "streaming_dir").mkdir()

    C.pipeline.enque(
        [
            C.fastq_parser.trim_reads(left="CCAGCC", right="GGTTCTGGC", tol=1),
            C.fastq_parser.translate(stop_readthrough=False),
            C.fastq_parser.save(where="pep", fmt="npy"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
        ]
    )

    streamer = C.data_loader.stream_from_fastq_dir(data_dir=data_dir)
    C.pipeline.stream(streamer=streamer, save_summary=False)
    
    lookup = Path(C.fastq_parser.dirs.logs)
    n_subs = sum(1 for p in lookup.iterdir() if p.is_dir())
    assert n_subs == 2
    return


def test_file_streaming(Walport_config, tmp_path, data_dir):
    import clibas as C

    C.initialize(Walport_config)

    C.fastq_parser.dirs.logs = str(tmp_path / "streaming_file")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "streaming_file")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "streaming_file")

    (tmp_path / "streaming_file").mkdir()

    C.pipeline.enque(
        [
            C.fastq_parser.translate(stop_readthrough=False),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
            C.fastq_parser.save(where="pep", fmt="npy"),
            C.fastq_parser.count_summary(where="pep", top_n=1000, fmt="csv"),
        ]
    )

    fname = Path(__file__).parent / "data" / "PADI4_r5_15_min_test.gz"
    streamer = C.data_loader.stream_from_gz_file(fname=fname, reads_per_chunk=500)
    C.pipeline.stream(streamer=streamer, save_summary=False)
    
    lookup = Path(C.fastq_parser.dirs.logs)
    n_subs = sum(1 for p in lookup.iterdir() if p.is_dir())
    assert n_subs == 4
    return

def test_streaming_memory_stable(tmp_path, data_dir):
    
    import gc
    import clibas as C

    C.initialize()

    C.fastq_parser.dirs.logs = str(tmp_path / "streaming_file")
    C.fastq_parser.dirs.parser_out = str(tmp_path / "streaming_file")
    C.fastq_parser.dirs.analysis_out = str(tmp_path / "streaming_file")

    (tmp_path / "streaming_file").mkdir()
    
    import tracemalloc
    tracemalloc.start()
    
    C.pipeline.enque([
        C.fastq_parser.translate(),
        C.fastq_parser.count_summary(where='pep', top_n=100, fmt='csv'),
    ])
    
    streamer = C.data_loader.stream_from_gz_dir(data_dir=data_dir)
    
    # measure memory before
    gc.collect()
    snapshot1 = tracemalloc.take_snapshot()
    
    # process streams
    C.pipeline.stream(streamer=streamer, save_summary=False)
    
    # measure after
    gc.collect()
    snapshot2 = tracemalloc.take_snapshot()
    
    # memory growth should be minimal (< 10MB for test data)
    stats = snapshot2.compare_to(snapshot1, 'lineno')
    total_growth = sum(stat.size_diff for stat in stats) / 1024 / 1024  # MB
    
    assert total_growth < 10, f"Memory leaked {total_growth:.1f}MB during streaming"
    tracemalloc.stop()























