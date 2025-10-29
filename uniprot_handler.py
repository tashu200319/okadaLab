"""
UniProtデータベースからの情報取得モジュール
"""

import re
import requests
from lxml import etree
import pandas as pd
from typing import List, Tuple, Optional


class UniprotData:
    """UniProtのXMLデータにアクセスし、情報を取得"""
    
    def __init__(self, uniprot_id: str):
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            self.xml = etree.fromstring(response.content)
            self.nsmap = self.xml.nsmap
            TF = self.xml.find('./', self.nsmap).text
            if TF != '\n  ':
                raise KeyError(f"Unexpected XML format: {TF}")
        except Exception as e:
            raise Exception(f"Failed to fetch UniProt data for {uniprot_id}: {e}")
    
    def get_id(self) -> List[str]:
        return [accession.text for accession in self.xml.findall('./entry/accession', self.nsmap)]
    
    def fasta(self) -> str:
        sequence = self.xml.find('./entry/sequence', self.nsmap)
        if sequence is not None:
            return sequence.text.replace('\n', '').replace(' ', '')
        return ""
    
    def get_fullname(self) -> str:
        fullname = self.xml.find('./entry/protein/*/fullName', self.nsmap)
        return fullname.text if fullname is not None else "No full name found"
    
    def get_organism(self) -> str:
        organism = self.xml.find('./entry/organism/name[@type="scientific"]', self.nsmap)
        return organism.text if organism is not None else "No organism found"
    
    def get_pdb_entries(self):
        pdb_entries = self.xml.findall(".//{http://uniprot.org/uniprot}dbReference[@type='PDB']", self.nsmap)
        return pdb_entries
    
    def getpdbdata(self, method) -> pd.DataFrame:
        if isinstance(method, str):
            methods = {m.strip() for m in re.split(r'[,\s]+', method) if m.strip()}
        else:
            methods = set(method)
        
        if not methods:
            methods = {"X-ray", "NMR", "EM"}
        
        pdbid = []
        data = []
        
        for dbReference in self.xml.findall('./entry/dbReference[@type="PDB"]', self.nsmap):
            x = []
            for propertys in dbReference:
                value = propertys.attrib["value"]
                x.append(value)
                if value == 'NMR':
                    x.append(None)
            
            if x and (x[0] in methods):
                pdbid.append(dbReference.attrib["id"])
                data.append(x)
        
        self.pdbdata = pd.DataFrame(data, index=pdbid, columns=['method', 'resolution', 'position']).T
        return self.pdbdata
    
    def pdblist(self, method="") -> List[str]:
        try:
            return self.pdbdata.columns.tolist()
        except AttributeError:
            return self.getpdbdata(method).columns.tolist()
    
    def position(self, pdbid: str) -> Tuple[int, int]:
        positiondata = self.pdbdata.at["position", pdbid].split(", ")
        
        if len(positiondata) == 1:
            _, posi = positiondata[0].split("=")
            beg, end = posi.split("-")
            beg = int(beg)
            end = int(end)
        else:
            beg = []
            end = []
            for position in positiondata:
                _, posi = position.split("=")
                align_beg, align_end = posi.split("-")
                beg.append(int(align_beg))
                end.append(int(align_end))
            beg = min(beg)
            end = max(end)
        
        return beg, end
