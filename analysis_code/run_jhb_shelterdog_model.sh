export STAN_NUM_THREADS=4
for i in {1..1}
	do
		time ./jhb-ordinal-hurdle-model sample \
		num_warmup=5000 num_samples=5000 save_warmup=0 \
		random seed=12345 id=$i data file=dog_data_for_cmdstan.R init=0 output file=jhb_dog_samples$i.csv &
	done
exit 0
