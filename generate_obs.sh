#!/bin/bash
mpirun -n 4 ./model -obs_gen T -obs_last 100 -obs_share 20
