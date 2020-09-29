from melissa_da_study import *
import os

preload_var="LD_PRELOAD"
preload_path=os.getcwd() + "/melissa-pdaf/install/lib/libpdaf_wrapper.so"
model_exe=os.getcwd() + "/model"

run_melissa_da_study(
        runner_cmd=model_exe,
        total_steps=100,
        ensemble_size=9,
        assimilator_type=ASSIMILATOR_PDAF,
        cluster=LocalCluster(),
        procs_server=4,
        procs_runner=2,
        n_runners=1,
        show_server_log = False,
        show_simulation_log = False,
        runner_timeout = 30,
        additional_server_env={
            preload_var: preload_path
            })

