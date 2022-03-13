import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
from pyteomics.usi import USI, proxi
from pyteomics.mass.unimod import Unimod
from ms2pip.single_prediction import SinglePrediction

def reshape_attributes(attributes):
    return {item["name"]: item["value"] for item in attributes}


def plot_unannotated_usi_spectrum(usi, peptide=None):
    proxi_response = proxi(usi)
    proxi_attributes = reshape_attributes(proxi_response["attributes"])

    # Process spectrum and plot
    spectrum = (
        sus.MsmsSpectrum(
            identifier=usi,
            precursor_mz=proxi_attributes["isolation window target m/z"],
            precursor_charge=int(proxi_attributes["charge state"]),
            mz=proxi_response["m/z array"],
            intensity=proxi_response["intensity array"],
            peptide=proxi_attributes["unmodified peptide sequence"],
            modifications=None
        )
        .filter_intensity(0.01, 50)
    )

    plt.figure(figsize=(12,6))
    plt.title(usi)
    sup.spectrum(spectrum)
    plt.show()


def plot_theoretical_spectrum(peptide): 
    ms2pip = SinglePrediction()
    mz, _, annotations = ms2pip.predict("VLHPLEGAVVIIFK", modifications="-", charge=2)
    #todo