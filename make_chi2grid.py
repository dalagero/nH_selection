import argparse
from dyb_analysis.fitter import fit
import prediction
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--save", required=True)
    parser.add_argument(
        "--set-up", choices=('rate', 'sin2', 'dm2', 'shape'), required=True
    )
    parser.add_argument("--sin22theta", type=float, required=True)
    parser.add_argument("--dm2ee", type=float, required=True)
    args = parser.parse_args()
    config_name=args.config
    theta13_values=0.5*np.arcsin(np.sqrt(np.array(args.sin22theta)))
    m2ee_values=np.array(args.dm2ee)
    constants=prediction.load_constants(config_name)
    starting_params=prediction.FitParams(0,0)
    if args.set_up == 'rate':
        frozen_params=fit.get_frozen_params(["all"],True,"pulled") #Rate-Only
        rate_only=True #Rate-only
    if args.set_up == 'sin2':
        frozen_params=fit.get_frozen_params(["all"],True,"free") #sin2-only scan
        rate_only=False #Rate+shape fits
    if args.set_up == 'dm2':
        frozen_params=fit.get_frozen_params(["all"],False,"frozen") #dm2-only scan
        rate_only=False #Rate+shape fits
    if args.set_up == 'shape':
        frozen_params=fit.get_frozen_params(["all"],True,"frozen") #Rate+shape 2-D
        rate_only=False #Rate+shape fits
    near_ads=None
    avg_near=True
    print("hereherehere")
    results=fit.grid(theta13_values,m2ee_values,constants,starting_params,frozen_params,near_ads,rate_only,avg_near)
    with np.load(args.save+".npz") as infile:
        theta13_array = infile['theta13_values']
        m2ee_array = infile['m2ee_values']
        chi2_values = infile['results']
    print(chi2_values)
    print(results)
    theta13_values = np.append(theta13_values,np.array(theta13_array))
    results = np.append(results,np.array(chi2_values))
    m2ee_values = np.append(m2ee_values,np.array(m2ee_array))
    np.savez(args.save,theta13_values=np.array(theta13_values),results=np.array(results),config=config_name,m2ee_values=np.array(m2ee_values),frozen_params=np.array(frozen_params))
    print(results)
