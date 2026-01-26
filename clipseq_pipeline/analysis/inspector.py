from bs4 import BeautifulSoup
from dataclasses import dataclass
import json
import requests
from typing import Callable

type InspectFn = Callable[[str, dict], dict]

def report_reader(filepath: str) -> None:
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data

def data_pubmed_search(organism: str, data: dict[str, dict[str, dict[str, float]]]) -> dict:
    proteins: dict = data["chi2_goodness_of_fit"]["results"]["mfe_lower_decile"]
    search_results: dict[str, dict] = {}
    protein_names: list[str] = []

    for name, stats in proteins.items():
        if name == "all_data":
            continue
        if stats["p-value"] > 0.05:
            continue

        search_result = search_pubmed(organism, name)
        combined_result = {**search_result, "stats": stats}

        search_results[name] = combined_result
        protein_names.append(name)

    sorted_results = dict(
        sorted(
            search_results.items(),
            key=lambda item: item[1].get("total_results", 0),
            reverse=True
        )
    )

    # Add the list of protein names to the dictionary
    return {
        "protein_names": protein_names,
        "results": sorted_results
    }

# TODO: refactor code.
def search_pubmed(organism: str, protein: str):
    # Step 1: Search for paper IDs
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "pubmed",
        "term": f"{organism} {protein.upper()} RNA motifs",
        "retmax": 10,  # top 10 papers
        "retmode": "json"
    }
    search_res = requests.get(search_url, params=search_params).json()
    id_list = search_res.get("esearchresult", {}).get("idlist", [])
    total_results = int(search_res.get("esearchresult", {}).get("count", 0))

    if not id_list:
        return {"total_results": total_results, "top_papers": []}

    # Step 2: Fetch paper details
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {
        "db": "pubmed",
        "id": ",".join(id_list),
        "retmode": "xml"
    }
    fetch_res = requests.get(fetch_url, params=fetch_params).text
    soup = BeautifulSoup(fetch_res, "lxml-xml")  # Use lxml-xml parser

    papers = []
    for article in soup.find_all("PubmedArticle"):
        title = article.ArticleTitle.text if article.ArticleTitle else "No title"
        pmid = article.PMID.text if article.PMID else ""
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}" if pmid else ""
        papers.append({"title": title, "url": url})

    return {"total_results": total_results, "top_papers": papers}

@dataclass(frozen=True)
class RnaInspection:
    name: str
    func: InspectFn

RNAMOTIFOLD_INSPECTION: list[RnaInspection] = [
    RnaInspection("paper_search_pubmed", data_pubmed_search)
]

def apply_inspection_report(
    organism: str,
    filepath: str,
    inspection_list: list[RnaInspection]
) -> None:
    """
    Run all inspection functions on the report and return combined results.
    """
    json_data = report_reader(f"{filepath}/report.json")
    results: dict = {}

    for inspection in inspection_list:
        func_result = inspection.func(organism, json_data)
        if not isinstance(func_result, dict):
            raise TypeError(f"Inspection function {inspection.name} must return a dict")
        results[inspection.name] = func_result
    if results:
        write_json(results, f"{filepath}/inspection.json")

def write_json(data: dict, filepath: str) -> None:
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)

def inspect_motifold_df(organism: str, filepath: str) -> None:
    return apply_inspection_report(organism, filepath, RNAMOTIFOLD_INSPECTION)


if __name__ == "__main__":
    apply_inspection_report("Caenorhabditis Elegans", "../../data/reports/caenorhabditis_elegans/rnamotifold", RNAMOTIFOLD_INSPECTION)