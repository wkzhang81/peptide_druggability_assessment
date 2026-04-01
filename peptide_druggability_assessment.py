#!/usr/bin/env python3
"""
Peptide Druggability Assessment Framework
All numerical outputs are standardized to two decimal places.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional
import warnings
import json
import os
from datetime import datetime

warnings.filterwarnings('ignore')

# RDKit Imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit not installed. Some molecular descriptor calculations will be limited.")
    RDKIT_AVAILABLE = False

# Plotting Configuration
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

class AminoAcidProperties:
    """Database of Amino Acid Physico-chemical Properties"""
    # Side chain pKa values (Empirical values)
    SIDE_CHAIN_PKA = {
        'D': 3.9, 'E': 4.3, 'H': 6.0, 'C': 8.3, 
        'Y': 10.1, 'K': 10.5, 'R': 12.5
    }
    
    # Kyte-Doolittle Hydrophobicity Index
    HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

class PeptideDruggabilityAssessor:
    """Core Class for Peptide Druggability Assessment"""
    
    def __init__(self, sequence: str, name: str = "Peptide_Target"):
        self.sequence = sequence.upper().strip()
        self.name = name
        self.aa_props = AminoAcidProperties()
        self.mol = None
        self.results = {}
        
        if RDKIT_AVAILABLE:
            self.mol = Chem.MolFromSequence(self.sequence)
            if self.mol:
                self.mol = Chem.AddHs(self.mol)

    def calculate_descriptors(self) -> Dict:
        """Calculate basic molecular descriptors (Two decimal places)"""
        if not RDKIT_AVAILABLE or not self.mol:
            return {"Error": "RDKit not available"}
        
        desc = {
            'MW': round(Descriptors.MolWt(self.mol), 2),
            'LogP': round(Descriptors.MolLogP(self.mol), 2),
            'TPSA': round(Descriptors.TPSA(self.mol), 2),
            'HBD': int(Descriptors.NumHDonors(self.mol)),
            'HBA': int(Descriptors.NumHAcceptors(self.mol)),
            'RotatableBonds': int(Descriptors.NumRotatableBonds(self.mol))
        }
        self.results['descriptors'] = desc
        return desc

    def assess_solubility(self, window_size: int = 5) -> Dict:
        """Solubility Assessment (CamSol-inspired sliding window)"""
        scores = []
        for i in range(len(self.sequence) - window_size + 1):
            window = self.sequence[i : i + window_size]
            h_sum = sum(self.aa_props.HYDROPHOBICITY.get(aa, 0) for aa in window)
            scores.append(h_sum / window_size)
        
        avg_hydro = np.mean(scores) if scores else 0
        camsol_like_score = 1 / (1 + np.exp(avg_hydro)) # Sigmoid transformation
        
        res = {
            'camsol_score': round(float(camsol_like_score), 2),
            'max_local_hydrophobicity': round(float(max(scores) if scores else 0), 2),
            'solubility_rank': "High" if camsol_like_score > 0.6 else "Medium" if camsol_like_score > 0.3 else "Low"
        }
        self.results['solubility'] = res
        return res

    def assess_charge_pka(self, ph_range: List[float] = None) -> Dict:
        """Calculate Net Charge via Henderson-Hasselbalch Equation"""
        if ph_range is None:
            ph_range = [2.0, 4.0, 7.0, 7.4, 10.0, 12.0]
        
        profile = {}
        for ph in ph_range:
            charge = 1.0 / (1.0 + 10**(ph - 9.6)) # N-term
            charge -= 1.0 / (1.0 + 10**(2.34 - ph)) # C-term
            
            for aa, pka in self.aa_props.SIDE_CHAIN_PKA.items():
                count = self.sequence.count(aa)
                if aa in 'DECY': 
                    charge -= count * (1.0 / (1.0 + 10**(pka - ph)))
                else: 
                    charge += count * (1.0 / (1.0 + 10**(ph - pka)))
            profile[round(ph, 2)] = round(charge, 2)
        
        self.results['charge_profile'] = profile
        return profile

    def assess_permeability(self) -> Dict:
        """Membrane Permeability Prediction"""
        desc = self.results.get('descriptors', self.calculate_descriptors())
        tpsa = desc.get('TPSA', 300)
        hbd = desc.get('HBD', 10)
        
        # Heuristic permeability score based on TPSA and HBD
        perm_score = (1 / (1 + (tpsa/150)**2)) * (1 / (1 + (hbd/8)**2))
        
        res = {
            'permeability_score': round(float(perm_score), 2),
            'penetration_potential': "Potential" if perm_score > 0.4 else "Weak"
        }
        self.results['permeability'] = res
        return res

    def assess_stability(self) -> Dict:
        """Metabolic Stability Prediction (Protease Cleavage Risk)"""
        cleavage_sites = 0
        for i in range(len(self.sequence) - 1):
            if self.sequence[i] in 'RK' and self.sequence[i+1] != 'P':
                cleavage_sites += 1
        
        stability_score = max(0.1, 1.0 - (cleavage_sites * 0.2))
        
        res = {
            'stability_score': round(float(stability_score), 2),
            'cleavage_sites_count': int(cleavage_sites),
            'risk_level': "High" if cleavage_sites > 3 else "Medium" if cleavage_sites > 1 else "Low"
        }
        self.results['stability'] = res
        return res

    def get_overall_score(self) -> float:
        """Calculate Final Druggability Score (Weighted Average)"""
        sol = self.results.get('solubility', {}).get('camsol_score', 0.5)
        perm = self.results.get('permeability', {}).get('permeability_score', 0.5)
        stab = self.results.get('stability', {}).get('stability_score', 0.5)
        
        # Weights: Stability (40%), Solubility (30%), Permeability (30%)
        overall = (sol * 0.3) + (perm * 0.3) + (stab * 0.4)
        final_score = round(float(overall), 2)
        self.results['overall_score'] = final_score
        return final_score

    def generate_report(self):
        """Generate Markdown Assessment Report"""
        score = self.results.get('overall_score', 0.0)
        report = f"# Peptide Druggability Assessment Report: {self.name}\n"
        report += f"- **Sequence**: {self.sequence}\n"
        report += f"- **Length**: {len(self.sequence)} AA\n"
        report += f"- **Overall Druggability Score**: {score:.2f}\n\n"
        
        report += "## Detailed Results\n"
        for module, data in self.results.items():
            if isinstance(data, dict):
                report += f"### {module.replace('_', ' ').capitalize()}\n"
                for k, v in data.items():
                    if isinstance(v, (float, np.float64)):
                        report += f"- {k}: {v:.2f}\n"
                    else:
                        report += f"- {k}: {v}\n"
                report += "\n"
        
        filename = f"{self.name}_report.md"
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"Report generated: {filename}")

    def plot_charge_profile(self):
        """Visualize pH-Charge Profile"""
        profile = self.results.get('charge_profile', {})
        if not profile: return
        
        plt.figure(figsize=(8, 5))
        x = list(profile.keys())
        y = list(profile.values())
        plt.plot(x, y, marker='o', linestyle='-', color='b')
        plt.axhline(0, color='red', linestyle='--')
        
        # Annotate points
        for xi, yi in zip(x, y):
            plt.text(xi, yi, f"{yi:.2f}", fontsize=9, verticalalignment='bottom')
            
        plt.title(f"pH-Net Charge Profile: {self.name}")
        plt.xlabel("pH")
        plt.ylabel("Net Charge")
        plt.grid(True, alpha=0.3)
        plt.savefig(f"{self.name}_charge_plot.png")
        plt.close()

def batch_process_to_csv(sequences: List[str], output_file: str = "peptide_results.csv"):
    """
    Processes a list of sequences and saves all assessment data into a CSV.
    """
    all_data = []

    for i, seq in enumerate(sequences):
        name = f"Peptide_{i+1}"
        print(f"Processing {name}: {seq}")
        
        # 1. Initialize and Run all assessments
        assessor = PeptideDruggabilityAssessor(seq, name=name)
        assessor.calculate_descriptors()
        assessor.assess_solubility()
        assessor.assess_charge_pka()
        assessor.assess_permeability()
        assessor.assess_stability()
        assessor.get_overall_score()
        
        # 2. Flatten the results dictionary for CSV format
        # We create a single-level dictionary for this row
        row = {
            "ID": name,
            "Sequence": seq,
            "Length": len(seq),
            "Overall_Score": assessor.results.get('overall_score')
        }
        
        # Add descriptors (MW, LogP, TPSA, etc.)
        row.update(assessor.results.get('descriptors', {}))
        
        # Add solubility metrics
        row.update(assessor.results.get('solubility', {}))
        
        # Add permeability and stability
        row.update(assessor.results.get('permeability', {}))
        row.update(assessor.results.get('stability', {}))
        
        # Add specific Net Charge at physiological pH (7.4)
        charge_profile = assessor.results.get('charge_profile', {})
        row['Net_Charge_pH7.4'] = charge_profile.get(7.4, "N/A")
        
        all_data.append(row)

    # 3. Create DataFrame and Export
    df = pd.DataFrame(all_data)
    
    # Reorder columns to make it readable (optional)
    cols = ['ID', 'Sequence', 'Overall_Score', 'MW', 'LogP', 'TPSA', 'solubility_rank', 'risk_level']
    # Add any remaining columns that weren't in the explicit list
    remaining_cols = [c for c in df.columns if c not in cols]
    df = df[cols + remaining_cols]
    
    df.to_csv(output_file, index=False)
    print(f"\nSuccessfully exported {len(sequences)} peptides to {output_file}")
