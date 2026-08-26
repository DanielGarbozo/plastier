import argparse, urllib.request, urllib.error, json, csv, time, sys

def fetch_json(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    return json.loads(urllib.request.urlopen(req).read().decode('utf-8'))

def get_ncbi_data(acc):
    try:
        uid = fetch_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={acc}&retmode=json")["esearchresult"]["idlist"][0]
        doc = fetch_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=json")["result"][uid]
        ftp = doc.get("ftppath_refseq") or doc.get("ftppath_genbank")
        if not ftp: return None, []

        rpt = urllib.request.urlopen(urllib.request.Request(f"{ftp}/{ftp.split('/')[-1]}_assembly_report.txt".replace("ftp://", "https://"))).read().decode('utf-8')
        plasmids, tot_bp, chr_bp, strain = [], 0, 0, doc.get("strain", "")
        
        for line in rpt.splitlines():
            if "strain=" in line: strain = line.split("strain=")[-1].strip()
            elif line.startswith("# Isolate:") and not strain: strain = line.split(":")[-1].strip()
            if line.startswith("#") or not line.strip(): continue
            p = line.split('\t')
            bp = int(p[8]) if p[8].isdigit() else 0
            tot_bp += bp
            if p[3].lower() == "chromosome": chr_bp += bp
            elif p[3].lower() == "plasmid": plasmids.append({"Plasmid_Accession": p[6] if p[6] != "na" else p[4], "Plasmid_Size_bp": bp})

        meta = {"Genome_Accession": acc, "Source_Link": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/", "Organism": doc.get("organism", ""), "Strain": strain, 
                "Assembly_Level": doc.get("assemblystatus", ""), "Total_Genome_Size_bp": tot_bp, 
                "Chromosome_Size_bp": chr_bp, "Plasmid_Count": len(plasmids)}
        return meta, plasmids
    except: 
        return None, []

def check_plsdb(acc):
    try:
        data = fetch_json(f"https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api/summary?NUCCORE_ACC={acc}")
        amr = data.get("Sequence_annotations", {}).get("AMR", [])
        g = list(set([r.get("gene_symbol") for r in amr if r.get("gene_symbol")]))
        c = list(set([r.get("drug_class") for r in amr if r.get("drug_class")]))
        return {"PLSDB_Found": "Yes", "AMR_Genes": "; ".join(g) or "None", "AMR_Classes": "; ".join(c) or "None"}
    except urllib.error.HTTPError as e: 
        return {"PLSDB_Found": "No" if e.code == 404 else f"Error {e.code}", "AMR_Genes": "N/A", "AMR_Classes": "N/A"}
    except: 
        return {"PLSDB_Found": "Error", "AMR_Genes": "Error", "AMR_Classes": "Error"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="table_summary.csv")
    args = parser.parse_args()

    with open(args.input) as f: accs = list(dict.fromkeys(f.read().split()))
    rows = []
    
    for acc in accs:
        print(f"Processing {acc}...")
        meta, plasmids = get_ncbi_data(acc)
        if not meta: continue
        
        if not plasmids:
            rows.append({**meta, "Plasmid_Accession": "None", "Plasmid_Size_bp": 0, "PLSDB_Found": "N/A", "AMR_Genes": "N/A", "AMR_Classes": "N/A"})
        
        for p in plasmids:
            amr = check_plsdb(p["Plasmid_Accession"])
            if amr["PLSDB_Found"] == "No" and p["Plasmid_Accession"].startswith("NZ_"): 
                amr = check_plsdb(p["Plasmid_Accession"][3:])
            rows.append({**meta, **p, **amr})
            time.sleep(1) # Be nice to PLSDB API
            
        time.sleep(0.5) # Be nice to NCBI API

    if rows:
        with open(args.output, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
        print(f"Saved {len(rows)} rows to {args.output}")
