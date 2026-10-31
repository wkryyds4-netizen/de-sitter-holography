"""Stage two: actual smeared de Sitter spectra, shared-clock two/four-point recovery."""
from pathlib import Path
import json, csv, hashlib, math
import numpy as np
import scipy
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_model import Spectrum, Smearing, matched_clocks

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/stage2"

def clock_description(name,clock):
    result={"family":name,"mean":clock.mean}
    if name != "exponential_energy_only":
        result.update(clock.record())
        result["time_tail_abs_gt_1"]=2*quad(lambda t:float(clock.time_density(t)),1,200,
                                           epsabs=1e-9,limit=1000)[0]
        result["time_tail_note"]="Numerical integral to |t|=200; residual checked in validation."
    else:
        result["time_rms"]=None
        result["time_note"]="Sharp positive-energy endpoint; not a matched finite-RMS clock."
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    records=[]
    spectra=[]
    for family in ("gaussian","cosine"):
        for sigma in (.1,.2,.3,.5):
            spectrum=Spectrum(Smearing(family,sigma),order=20,periods=50)
            spectra.append(spectrum.record())
            tasks=[(mean,name,clock) for mean in (5.,10.,20.)
                   for name,clock in matched_clocks(mean,.5).items()]
            errors=spectrum.errors_many([task[2] for task in tasks])
            for (mean,name,clock),error in zip(tasks,errors):
                row={"field_family":family,"field_sigma":sigma,"clock_mean":mean,
                     "clock_family":name,"clock":clock_description(name,clock),
                     "spectrum":spectrum.record(),**error}
                records.append(row)
            print(f"Completed {family}, sigma={sigma}",flush=True)
    data={"scope":"Normalized absolute error, not pointwise relative OTOC error. Fixed dS4 free field, no recoil or finite G.",
          "normalization":"A=phi_f/sqrt(Z), C2_bulk(0)=1, C4_bulk(0,0,0,0)=3",
          "spectra":spectra,"cases":records,
          "versions":{"numpy":np.__version__,"scipy":scipy.__version__},
          "code_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in [Path(__file__),ROOT/"code/field_model.py",ROOT/"code/clock_model.py"]}}
    (OUT/"study.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    fields=["field_family","field_sigma","clock_mean","clock_family","sup_C2_absolute",
            "sup_C4_absolute_over_3","dressed_fourth_cumulant","spectral_truncation_C4_over_3_bound"]
    with (OUT/"errors.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        writer.writerows({k:r[k] for k in fields} for r in records)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.5),constrained_layout=True)
    for ax,family in zip(axes,("gaussian","cosine")):
        for name,color in (("compact","#126d82"),("gamma5","#7350a2"),("exponential_energy_only","#b24935")):
            selected=[r for r in records if r["field_family"]==family and r["clock_family"]==name and r["clock_mean"]==10]
            x=[r["field_sigma"] for r in selected]
            ax.semilogy(x,[r["sup_C2_absolute"] for r in selected],color=color,ls="--",marker="o")
            ax.semilogy(x,[r["sup_C4_absolute_over_3"] for r in selected],color=color,ls="-",marker="s",label=name)
        ax.set(xlabel="Field smearing RMS in modular time",ylabel="Worst absolute error / zero-time bulk scale",
               title=f"{family}: solid=four-point, dashed=two-point")
        ax.legend(fontsize=8)
    fig.suptitle("Physical dS4 field spectra: two-point accuracy can hide larger four-point errors\nMean clock energy=10; compact and Gamma have time RMS=0.5; exponential matches energy only",fontsize=11)
    fig.savefig(OUT/"two-vs-four.png",dpi=180);fig.savefig(OUT/"two-vs-four.svg");plt.close(fig)

    spec=Spectrum(Smearing("gaussian",.3),order=20)
    clock=matched_clocks(10,.5)["compact"]
    times=np.linspace(0,3,61)
    curve=[]
    for t in times:
        c2=spec.c2(np.array(t))
        bulk=complex(c2*c2+1+abs(c2)**2)
        dressed=spec.c4((t,0,t,0),clock)
        curve.append({"t":float(t),"bulk_real":bulk.real,"bulk_imag":bulk.imag,
                      "clock_real":dressed.real,"clock_imag":dressed.imag,
                      "absolute_error_over_3":abs(dressed-bulk)/3})
    (OUT/"otoc-curve.json").write_text(json.dumps(curve,indent=2)+"\n",encoding="utf-8")
    bound=spec.errors(clock)["sup_C4_absolute_over_3"]
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    axes[0].plot(times,[r["bulk_real"] for r in curve],label="Bulk free-field Re C4")
    axes[0].plot(times,[r["clock_real"] for r in curve],label="Shared-clock Re C4")
    axes[0].set(xlabel="Modular time t",ylabel="Literal A(t) A(0) A(t) A(0)",title="A dressed Gaussian field need not have Gaussian moments")
    axes[0].legend(fontsize=8)
    axes[1].plot(times,[r["absolute_error_over_3"] for r in curve],color="#126d82")
    axes[1].axhline(bound,color="#b24935",ls="--",label="All-real-time supremum (attained at t=0)")
    axes[1].set(xlabel="Modular time t",ylabel="Absolute four-point error / 3",title="No time-growth assumption was imposed")
    axes[1].legend(fontsize=8)
    fig.savefig(OUT/"shared-clock-otoc.png",dpi=180);fig.savefig(OUT/"shared-clock-otoc.svg");plt.close(fig)
    print(json.dumps({"cases":len(records),"reference":next(r for r in records if r["field_family"]=="gaussian"
                        and r["field_sigma"]==.3 and r["clock_mean"]==10 and r["clock_family"]=="compact")},indent=2))

if __name__=="__main__":
    main()

