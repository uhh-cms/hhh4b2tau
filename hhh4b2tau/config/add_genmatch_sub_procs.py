import order as od
import law

logger = law.logger.get_logger(__name__)


def add_genmatch_subprocesses(cfg):
    logger.info("adding sub processes for genmatching")
    for proc in cfg.processes.values():
        if proc.has_tag("hhh"):
            try:
                proc.add_process(
                    name=f"{proc.name}_nomatch",
                    id=proc.id + 100,
                    label=f"{proc.label} unmatched",
                )

                # super process if h1 is matched
                h1_proc = proc.add_process(
                    name=f"{proc.name}_match_h1",
                    id=proc.id + 200,
                    label=f"{proc.label} match $H_{{1}}$",
                )
                proc_h1_no_h2 = h1_proc.add_process(
                    name=f"{proc.name}_match_h1_no_h2",
                    id=proc.id + 300,
                    label=f"{proc.label} match $H_{{1}}, !H_{{2}}$",
                )
                proc_h1_and_h2 = h1_proc.add_process(
                    name=f"{proc.name}_match_h1_and_h2",
                    id=proc.id + 400,
                    label=f"{proc.label} match $H_{{1}} + H_{{2}}$",
                )
                h2_proc = proc.add_process(
                    name=f"{proc.name}_match_h2",
                    id=proc.id + 500,
                    label=f"{proc.label} match $H_{{2}}$",
                )
                proc_h2_no_h1 = h2_proc.add_process(
                    name=f"{proc.name}_match_h2_no_h1",
                    id=proc.id + 600,
                    label=f"{proc.label} match $H_{{2}}, !H_{{1}}$",
                )
                h2_proc.add_process(proc_h1_and_h2)
            except od.unique.DuplicateNameException as e:
                from IPython import embed
                embed(header=f"raised exception {e}")
                raise e

