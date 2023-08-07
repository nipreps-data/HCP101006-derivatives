## Fieldmap coefficients file

Results here were calculated with topup and the LR/RL epi blips under the HCP1011006/fmap of the test dataset, using the following script and SDCFlows 2.5.2.dev47:

```Python
from pathlib import Path
import json
import numpy as np
from nipype import config
from nipype.pipeline import engine as pe
from sdcflows import fieldmaps as sfm
from sdcflows import transform as st
from sdcflows.workflows.apply import correction as swac
from sdcflows.workflows.apply.registration import init_coeff2epi_wf
from mriqc.workflows.anatomical.base import synthstrip_wf

pe0 = "LR"
pe1 = pe0[::-1]

config.set("execution", "remove_unnecessary_outputs", "false")

datapath = Path("/data/datasets/sdcflows-tests/HCP101006/sub-101006")

metadata_pe0 = json.loads((datapath / "fmap" / f"sub-101006_dir-{pe0}_epi.json").read_text())
metadata_pe1 = json.loads((datapath / "fmap" / f"sub-101006_dir-{pe1}_epi.json").read_text())
epi_blips = [
    str(datapath / "fmap" / f"sub-101006_dir-{pe0}_epi.nii.gz"),
    str(datapath / "fmap" / f"sub-101006_dir-{pe1}_epi.nii.gz"),
]

estimator = sfm.FieldmapEstimation(sources=epi_blips, bids_id=f"hcp_pepolar{pe0}")

# Enable this and disable the inputs for phasediff estimation
# estimator = sfm.FieldmapEstimation(sources=[f"data/sub-pilot_ses-{session}_phasediff.nii.gz"], bids_id="phasediff0")

step1 = estimator.get_workflow(omp_nthreads=12, debug=False)
step1.inputs.inputnode.metadata = [metadata_pe0, metadata_pe1]
step1.inputs.inputnode.in_data = epi_blips

bex_sbref = synthstrip_wf(name="bex_sbref_wf", omp_nthreads=12)
bex_sbref.inputs.inputnode.in_files = str(
    datapath / "func" / f"sub-101006_task-rest_dir-{pe0}_sbref.nii.gz"
)
align_wf = init_coeff2epi_wf(omp_nthreads=12)

step2 = swac.init_unwarp_wf(omp_nthreads=12)
step2.inputs.inputnode.distorted = str(
    datapath / "func" / f"sub-101006_task-rest_dir-{pe0}_sbref.nii.gz"
)
step2.inputs.inputnode.metadata = json.loads(
    (datapath / "func" / f"sub-101006_task-rest_dir-{pe0}_sbref.json").read_text()
)

skull_stripping = synthstrip_wf(omp_nthreads=12)

wf = pe.Workflow(name=f"hcp_testsdc{pe0}")
wf.base_dir = "work/"
wf.connect([
    (step1, step2, [
        ("outputnode.fmap_coeff", "inputnode.fmap_coeff"),
    ]),
    (step1, align_wf, [
        ("outputnode.fmap_ref", "inputnode.fmap_ref"),
        ("outputnode.fmap_mask", "inputnode.fmap_mask"),
    ]),
    (bex_sbref, align_wf, [
        ("outputnode.out_corrected", "inputnode.target_ref"),
        ("outputnode.out_mask", "inputnode.target_mask"),
    ]),
    (align_wf, step2, [
        ("outputnode.target2fmap_xfm", "inputnode.data2fmap_xfm")
    ]),
    (step2, skull_stripping, [
        ("outputnode.corrected_ref", "inputnode.in_files"),
    ]),
])

wf.run()
```

File mappings:

* `sub-101006_coeff-1_desc-topup_fieldmap.nii.gz` <- `<work_dir>/hcp_testsdcLR/wf_hcp_pepolarLR/fix_coeff/sub-101006_dir-RL_epi_average_merged_sliced_volreg_base_fieldcoef_fixed.nii.gz`.
* `sub-101006_desc-pepolar_epiref.nii.gz` <- `<work_dir>/hcp_testsdcLR/wf_hcp_pepolarLR/brainextraction_wf/n4/clipped_corrected.nii.gz`.
* `sub-101006_desc-pepolar_mask.nii.gz` <- `<work_dir>/hcp_testsdcLR/wf_hcp_pepolarLR/brainextraction_wf/masker/clipped_mask.nii.gz`.
* `sub-101006_task-rest_dir-LR_desc-preproc_sbref.nii.gz` <- `<work_dir>/hcp_testsdcLR/bex_sbref_wf/post_n4/clipped_corrected.nii.gz`.
* `sub-101006_task-rest_dir-LR_desc-sbref_mask.nii.gz` <- `<work_dir>/hcp_testsdcLR/bex_sbref_wf/synthstrip/clipped_corrected_desc-brain_mask.nii.gz`.
* `sub-101006_from-sbrefLR_to-fieldmapref_mode-image_xfm.mat` <-`<work_dir>/hcp_testsdcLR/coeff2epi_wf/coregister/transform0GenericAffine.mat`