"""
UniProtデータベースからの情報取得モジュール (XML/JSONフォールバック版)
"""

import re
import json
import requests
from lxml import etree
import pandas as pd
from typing import List, Tuple, Optional


class UniprotData:
    """UniProtのデータにアクセスし、情報を取得（XML優先、失敗時JSON）"""
    
    def __init__(self, uniprot_id: str):
        self.uniprot_id = uniprot_id
        self.data_format = None  # 'xml' or 'json'
        self.xml = None
        self.json_data = None
        self.nsmap = None
        
        # まずXML形式を試す
        if self._try_xml(uniprot_id):
            self.data_format = 'xml'
        # XML失敗時はJSON形式を試す
        elif self._try_json(uniprot_id):
            self.data_format = 'json'
        else:
            raise Exception(f"Failed to fetch UniProt data for {uniprot_id} in both XML and JSON formats")
    
    def _try_xml(self, uniprot_id: str) -> bool:
        """XML形式でデータ取得を試行"""
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            self.xml = etree.fromstring(response.content)
            self.nsmap = self.xml.nsmap
            
            # XMLフォーマットのチェック
            TF = self.xml.find('./', self.nsmap)
            if TF is not None and TF.text and TF.text != '\n  ':
                # 著作権メッセージなどが返ってきた場合
                if 'Copyright' in TF.text or 'Creative Commons' in TF.text:
                    return False
            
            # entry要素が存在するかチェック
            entry = self.xml.find('./entry', self.nsmap)
            if entry is None:
                return False
            
            return True
            
        except Exception as e:
            # XMLでの取得失敗
            return False
    
    def _try_json(self, uniprot_id: str) -> bool:
        """JSON形式でデータ取得を試行"""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=json"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            self.json_data = response.json()
            return True
        except Exception as e:
            return False
    
    def get_id(self) -> List[str]:
        """UniProt IDのリストを取得"""
        if self.data_format == 'xml':
            return [accession.text for accession in self.xml.findall('./entry/accession', self.nsmap)]
        else:  # json
            ids = [self.json_data.get('primaryAccession', '')]
            ids.extend(self.json_data.get('secondaryAccessions', []))
            return [id for id in ids if id]
    
    def fasta(self) -> str:
        """FASTA配列を取得"""
        if self.data_format == 'xml':
            sequence = self.xml.find('./entry/sequence', self.nsmap)
            if sequence is not None:
                return sequence.text.replace('\n', '').replace(' ', '')
            return ""
        else:  # json
            seq_data = self.json_data.get('sequence', {})
            return seq_data.get('value', '').replace('\n', '').replace(' ', '')
    
    def get_fullname(self) -> str:
        """タンパク質のフルネームを取得"""
        if self.data_format == 'xml':
            fullname = self.xml.find('./entry/protein/*/fullName', self.nsmap)
            return fullname.text if fullname is not None else "No full name found"
        else:  # json
            protein_desc = self.json_data.get('proteinDescription', {})
            recommended_name = protein_desc.get('recommendedName', {})
            if recommended_name:
                full_name = recommended_name.get('fullName', {})
                return full_name.get('value', 'No full name found')
            return "No full name found"
    
    def get_organism(self) -> str:
        """生物種を取得"""
        if self.data_format == 'xml':
            organism = self.xml.find('./entry/organism/name[@type="scientific"]', self.nsmap)
            return organism.text if organism is not None else "No organism found"
        else:  # json
            organism = self.json_data.get('organism', {})
            return organism.get('scientificName', 'No organism found')
    
    def get_pdb_entries(self):
        """PDBエントリを取得"""
        if self.data_format == 'xml':
            return self.xml.findall(".//{http://uniprot.org/uniprot}dbReference[@type='PDB']", self.nsmap)
        else:  # json
            # JSON形式の場合は別の形式で返す
            return self._get_pdb_entries_from_json()
    
    def _get_pdb_entries_from_json(self):
        """JSON形式からPDBエントリを抽出（内部用）"""
        pdb_entries = []
        cross_refs = self.json_data.get('uniProtKBCrossReferences', [])
        for ref in cross_refs:
            if ref.get('database') == 'PDB':
                pdb_entries.append(ref)
        return pdb_entries
    
    def _normalize_method_label(self, method_val: str) -> Optional[str]:
        """
        UniProt JSON の Method 文字列を
        'X-ray' / 'EM' / 'NMR' / None に正規化する。
        """
        if not method_val:
            return None

        s = method_val.lower()

        # X-ray 系
        if "x-ray" in s or "xray" in s:
            return "X-ray"

        # EM 系
        if "electron microscopy" in s or "electron cryo-microscopy" in s or "cryo-em" in s:
            return "EM"

        # NMR 系
        if "nmr" in s:
            return "NMR"

        # 対応外の Method は無視
        return None
    
    
    def getpdbdata(self, method) -> pd.DataFrame:
        """PDBデータを取得"""
        if isinstance(method, str):
            methods = {m.strip() for m in re.split(r'[,\s]+', method) if m.strip()}
        else:
            methods = set(method)
        
        if not methods:
            methods = {"X-ray", "NMR", "EM"}
        
        pdbid = []
        data = []
        
        if self.data_format == 'xml':
            # XML形式の処理（既存のロジック）
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
        
        else:  # json
            # JSON形式の処理
            cross_refs = self.json_data.get('uniProtKBCrossReferences', [])
            for ref in cross_refs:
                if ref.get('database') != 'PDB':
                    continue

                pdb_id = ref.get('id')
                properties = ref.get('properties', [])

                raw_method = None
                resolution_val = None
                chains_val = None

                # プロパティから必要な情報を抽出
                for prop in properties:
                    key = prop.get('key')
                    value = prop.get('value')

                    if key == 'Method':
                        raw_method = value
                    elif key == 'Resolution':
                        resolution_val = value
                    elif key == 'Chains':
                        chains_val = value

                # Method を正規化（X-ray / EM / NMR / None）
                norm_method = self._normalize_method_label(raw_method)
                if not norm_method:
                    # 対応外の Method はスキップ
                    continue

                # 呼び出し側から method 指定がある場合だけフィルタ
                # 例: methods = {"X-ray", "EM"}
                if methods and norm_method not in methods:
                    continue

                x = [norm_method, resolution_val, chains_val]
                if norm_method == 'NMR':
                    # NMR のときは resolution を None にするという元の仕様を踏襲
                    x[1] = None

                pdbid.append(pdb_id)
                data.append(x)

        
        self.pdbdata = pd.DataFrame(data, index=pdbid, columns=['method', 'resolution', 'position']).T
        return self.pdbdata
    
    def pdblist(self, method="") -> List[str]:
        """PDB IDのリストを取得"""
        try:
            return self.pdbdata.columns.tolist()
        except AttributeError:
            return self.getpdbdata(method).columns.tolist()
    
    def position(self, pdbid: str) -> Tuple[Optional[int], Optional[int]]:
        """
        PDB IDの位置情報を取得
        
        Parameters
        ----------
        pdbid : str
            PDB ID
        
        Returns
        -------
        tuple
            (beg, end) または (None, None)
        """
        try:
            # pdbdataから位置情報を取得
            position_value = self.pdbdata.at["position", pdbid]
            
            # Noneまたは空チェック
            if position_value is None or pd.isna(position_value):
                print(f"  Warning: No position data for PDB {pdbid}")
                return None, None
            
            # 文字列でない場合
            if not isinstance(position_value, str):
                print(f"  Warning: Invalid position data type for PDB {pdbid}: {type(position_value)}")
                return None, None
            
            # 位置情報をパース
            positiondata = position_value.split(", ")
            
            if len(positiondata) == 1:
                # 単一チェーンの場合
                if "=" not in positiondata[0]:
                    print(f"  Warning: Invalid position format for PDB {pdbid}: {positiondata[0]}")
                    return None, None
                
                _, posi = positiondata[0].split("=")
                if "-" not in posi:
                    print(f"  Warning: Invalid position range for PDB {pdbid}: {posi}")
                    return None, None
                
                beg, end = posi.split("-")
                beg = int(beg)
                end = int(end)
            else:
                # 複数チェーンの場合
                beg = []
                end = []
                for position in positiondata:
                    if "=" not in position:
                        print(f"  Warning: Skipping invalid position format: {position}")
                        continue
                    
                    _, posi = position.split("=")
                    if "-" not in posi:
                        print(f"  Warning: Skipping invalid position range: {posi}")
                        continue
                    
                    align_beg, align_end = posi.split("-")
                    beg.append(int(align_beg))
                    end.append(int(align_end))
                
                if not beg or not end:
                    print(f"  Warning: No valid position data found for PDB {pdbid}")
                    return None, None
                
                beg = min(beg)
                end = max(end)
            
            return beg, end
            
        except KeyError:
            print(f"  Warning: PDB {pdbid} not found in pdbdata")
            return None, None
        except (ValueError, AttributeError) as e:
            print(f"  Warning: Error parsing position for PDB {pdbid}: {e}")
            return None, None
        except Exception as e:
            print(f"  Error: Unexpected error getting position for PDB {pdbid}: {e}")
            return None, None