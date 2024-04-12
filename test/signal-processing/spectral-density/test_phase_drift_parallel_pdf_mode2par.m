% 2022-11-09 12:41:41.357980800 +0100

	Scy = 1;

	sy = phase_drift_parallel_pdf_mode2par(Scy)

	Scy_ = phase_drift_parallel_pdf(0,sy)

	fail = abs(Scy_ - Scy) > sqrt(eps)

