bsub -I -b -q q_sw_expr -host_stack 1024 -share_size 4096 -n 10 -cgsp 64 ./a.out
