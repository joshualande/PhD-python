from os.path import join, basename
from textwrap import dedent
from itertools import product

from lande.utilities.jobtools import JobBuilder

class PipelineBuilder(JobBuilder):
    def build_followup(self, code, hypotheses, followups, followup_filenames=None, extra=''):
        if followup_filenames is None: followup_filenames = followups
        for perm in product(*self.values):
            subdir = self.build_subdir(perm)
            args = self.build_args(perm)

            for hypothesis in hypotheses:
                for followup,followup_filename in zip(followups,followup_filenames):
                    pwn = basename(subdir)
                    run = join(subdir,'followup_%s_%s_%s.sh' % (pwn,followup_filename,hypothesis))
                    open(run,'w').write("python %s %s %s --hypothesis=%s --followup=%s" % (code,args,extra,hypothesis,followup))

            followup_all = join(self.savedir,'followup_all.sh')
            open(followup_all,'w').write(dedent("""
            #!/usr/bin/env bash
            submit_list=""
            for pwn in `find . -maxdepth 1 -type d|xargs -i basename {}`; do
                for hypothesis in at_pulsar point extended; do
                    if [ -e $pwn/roi_${hypothesis}_${pwn}.dat -a -e $pwn/results_${pwn}_pointlike_${hypothesis}.yaml ]; then
                        submit_list="$pwn/followup_${pwn}_*_${hypothesis}.sh $submit_list"
                    fi
                done
            done
            submit_all $submit_list $@
            """))

