import os
import sys
import argparse
import multiprocessing

from math import ceil

from readfq import readfq # cheers heng

def variant_worker(write_q, ref_seq, msa, window_i, start_record_i, end_record_i):
    # Open the MSA, iterate over each sequence and walk the genome to find
    # diagreements with the loaded reference
    # NOTE This particular MSA does not handle insertions
    record_i = -1
    first = True
    with open(msa, 'r') as all_fh:
        for name, seq, qual in readfq(all_fh):
            record_i += 1
            if record_i < start_record_i:
                continue
            if record_i > end_record_i:
                break

            if first:
                sys.stderr.write("[STAT] Worker %d started on record %d\n" % (window_i, record_i))
                first = False

            central_sample_id = name

            query_on_ref_pos = 0
            current_deletion_len = 0

            curr_lines = []
            for qbase in seq:
                if qbase == '-':
                    # Extend the length of the current deletion
                    current_deletion_len += 1
                else:
                    if current_deletion_len > 0:
                        # We've come to the end of a deletion, output it
                        curr_lines.append(','.join([
                            central_sample_id,
                            #str( "%d-%d" % ((query_on_ref_pos-current_deletion_len)+1, query_on_ref_pos) ),
                            str((query_on_ref_pos-current_deletion_len)+1),
                            "",
                            "%dD" % current_deletion_len,
                            "1",
                        ]))
                        current_deletion_len = 0

                # Now deletions are handled, check for single nucleotide variants
                # NOTE This includes missing data such as N
                # NOTE This algorithm does not consider INS against ref
                if qbase != ref_seq[query_on_ref_pos]:
                    if current_deletion_len == 0:
                        # SNV detected and we aren't in an active DEL
                        curr_lines.append(','.join([
                            central_sample_id,
                            str(query_on_ref_pos+1),
                            ref_seq[query_on_ref_pos],
                            qbase,
                            "0",
                        ]))

                # Advance pointer (this is overkill here but a useful starting point
                # for a future algo walking the ref for insertions)
                query_on_ref_pos += 1

            if current_deletion_len > 0:
                # Output the last deletion, if there is one 
                # (this is almost always going to be garbage but we include it for completeness)
                curr_lines.append(','.join([
                    central_sample_id,
                    #str( "%d-%d" % ((query_on_ref_pos-current_deletion_len)+1, query_on_ref_pos) ),
                    str((query_on_ref_pos-current_deletion_len)+1),
                    "",
                    "%dD" % current_deletion_len,
                    "1",
                ]))

            # Push curr lines to writer
            write_q.put( '\n'.join(curr_lines) + '\n' )

        # Break out, send sentinel to queue
        sys.stderr.write("[DONE] Worker %d finished at next record %d\n" % (window_i, record_i))
        write_q.put(None)

def write_worker(write_q, n):
    dead_workers = 0
    sys.stdout.write(','.join([
        "COG-ID",
        "Position",
        "Reference_Base",
        "Alternate_Base",
        "Is_Indel"
    ]) + '\n')
    while True:
        work = write_q.get()
        if not work:
            dead_workers += 1
            sys.stderr.write("[NOTE] Writer closed connection to worker. Waiting on %d more workers.\n" % (n - dead_workers))
        else:
            sys.stdout.write(work)

        if dead_workers == n:
            sys.stderr.write("[NOTE] No workers remaining. Closing writer.\n")
            return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--msa", required=True)
    parser.add_argument("-n", type=int, required=True, help="number of sequences to process, will be divided amongst threads")
    parser.add_argument("-t", "--threads", type=int, default=4)
    args = parser.parse_args()

    # Check files exist
    for fpt, fp in ("REF", args.ref), ("MSA", args.msa):
        if not os.path.isfile(fp):
            sys.stderr.write("[FAIL] Could not open %s %s.\n" % (fpt, fp))
            sys.exit(1)
        else:
            sys.stderr.write("[NOTE] %s: %s\n" % (fpt, fp))

    sys.stderr.write("[NOTE] NUM_SEQUENCES: %d\n" % args.n)
    sys.stderr.write("[NOTE] ASKLEPIAN_VARIANT_THREADS: %d\n" % args.threads)

    # Load the ref and assign it to ref_seq
    with open(args.ref) as canon_fh:
        for name, seq, qual in readfq(canon_fh):
            break
        if not name:
            sys.stderr.write("[FAIL] Could not read sequence from reference.\n")
            sys.exit(2)
        else:
            ref_seq = seq

    write_q = multiprocessing.Queue()
    processes = []

    writer_process = multiprocessing.Process(
        target=write_worker,
        args=(
            write_q,
            args.threads,
        ),
    )
    processes.append(writer_process)

    window_l = ceil(args.n / float(args.threads))
    for window_i, window_pos in enumerate(range(0, args.n, window_l)):
        start = window_pos
        end = window_pos + window_l - 1 # remove 1 as we dont use gte in worker
        if window_i == (args.threads - 1):
            end = args.n # in case we've managed to screw the last window and its too short, just set it to N

        sys.stderr.write("[WORK] Worker %d (%d, %d)\n" % (window_i, start, end))
        p = multiprocessing.Process(
            target=variant_worker,
            args=(
                write_q,
                ref_seq,
                args.msa,
                window_i,
                start,
                end,
            ),
        )
        processes.append(p)

    # Engage
    for p in processes:
        p.start()

    # Block
    for p in processes:
        p.join()

    sys.stderr.write("[DONE] All workers exited, bye!\n")


if __name__ == "__main__":
    main()
