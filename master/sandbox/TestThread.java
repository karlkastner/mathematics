	class TestThread
	{
		int MAX_NUM_THREADS=8;
		int tnum;
		Worker[] t = new Worker[MAX_NUM_THREADS];
		int n; //= 40000;
		int [] s = new int[8];
		int [] p = new int[8];

		int P[];

		public TestThread() {
				P = new int[2];
				P[0] = 0;
				P[1] = 0;
			};
	
		public void testThread(int n_) throws InterruptedException
		{
			n = n_;

			// query the number of requested threads
			String omp_n_thread_s = System.getenv("OMP_NUM_THREADS");
			int omp_n_thread;
			try {
				omp_n_thread = Integer.parseInt( omp_n_thread_s );
			} catch ( Exception e )
			{
				omp_n_thread = java.lang.Integer.MAX_VALUE;
			}

			// query the number of available CPUs
			int available_processors = Math.max(1, java.lang.Runtime.getRuntime().availableProcessors());

			// determine the number of threads
			tnum = Math.min(Math.min( MAX_NUM_THREADS, available_processors), omp_n_thread);

			// Why do the workers have to be created twice?
			for (int i=0; i<tnum; i++)
			{	
				t[i] = new Worker();
			}
			for (int i=0; i<tnum; i++)
			{
				t[i].start();
			}
			for (int i=0; i<tnum; i++)
			{
				t[i].join();
			}
	
			int ss=0;
			int pp=0;
			for (int i=0; i<tnum; i++)
			{
				ss = ss+s[i];
				pp = pp + p[i];
			}	
			System.out.println(ss + " " + pp);
		} // testThread()

		class Worker extends Thread
		{
			public synchronized void run()
			{
				for (int i=0; i<tnum; i++)
				{
					if (t[i] == this)
					{
						s[i] = 0;
						for (int j=0; j<n/tnum; j++)
						{
							//int p = 0;
							for (int k=0; k<1000000; k++)
							{
								//p[i]=p[i]+1;
								p[P[i]]=p[P[i]]+1;
							}
							s[i]=s[i]+1;
							//System.out.println(p); //getId() + " " + i + " " + tnum);
						}
						//System.out.println(getId() + " " + i + " " + tnum);
					}
				}
			} // run()
		} // class Worker
	} // TestThread

