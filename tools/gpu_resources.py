"""Conservative GPU availability checks; never preempt an existing process."""
import csv,io,subprocess
def snapshot():
    text=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True)
    processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader,nounits'],text=True)
    occupied={r[0].strip() for r in csv.reader(io.StringIO(processes)) if r}
    return [dict(gpu=int(r[0]),uuid=r[1].strip(),memory_MiB=int(r[2]),utilization=int(r[3]),has_compute_process=r[1].strip() in occupied) for r in csv.reader(io.StringIO(text)) if r]
def assert_gpu_free(gpu):
    r=next(r for r in snapshot() if r['gpu']==gpu)
    assert r['memory_MiB']<100 and not r['has_compute_process'] and r['utilization']==0,('GPU unavailable; do not interfere with existing work',r)
    return r
