import urllib.request
import urllib
# queryText = """

# <?xml version="1.0" encoding="UTF-8"?>

# <orgPdbQuery>
#     <version>head</version>
#     <queryType>org.pdb.query.simple.NoLigandQuery</queryType>
#     <description>Ligand Search : Has free ligands=no</description>
#     <queryId>2AC020D3</queryId>
#     <resultCount>36116</resultCount>
#     <runtimeStart>2019-01-16T17:04:56Z</runtimeStart>
#     <runtimeMilliseconds>1817</runtimeMilliseconds>
#     <haveLigands>no</haveLigands>

#     <version>head</version>
#     <queryType>org.pdb.query.simple.SecondaryStructureQuery</queryType>
#     <description>Secondary structure has:  0 or less Beta Sheets</description>
#     <queryId>79E5045E</queryId>
#     <resultCount>14916</resultCount>
#     <runtimeStart>2019-01-16T17:35:18Z</runtimeStart>
#     <runtimeMilliseconds>402</runtimeMilliseconds>
#     <polyStats.helixPercent.comparator>between</polyStats.helixPercent.comparator>
#     <polyStats.helixCount.comparator>between</polyStats.helixCount.comparator>
#     <polyStats.sheetPercent.comparator>between</polyStats.sheetPercent.comparator>
#     <polyStats.sheetCount.comparator>between</polyStats.sheetCount.comparator>
#     <polyStats.sheetCount.max>0</polyStats.sheetCount.max>
#   </orgPdbQuery>
# """
queryText = """
<orgPdbCompositeQuery version="1.0">
    <resultCount>36116</resultCount>
    <queryId>B12E5C49</queryId>
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.NoLigandQuery</queryType>
    <description>Ligand Search : Has free ligands=no</description>
    <queryId>D1BA680C</queryId>
    <resultCount>36116</resultCount>
    <runtimeStart>2019-01-16T17:38:27Z</runtimeStart>
    <runtimeMilliseconds>427</runtimeMilliseconds>
    <haveLigands>no</haveLigands>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.LargeStructureQuery</queryType>
    <description>Search for Large Structures : search option = omit large structures</description>
    <queryId>7A9CBD68</queryId>
    <resultCount>147263</resultCount>
    <runtimeStart>2019-01-16T17:42:56Z</runtimeStart>
    <runtimeMilliseconds>1593</runtimeMilliseconds>
    <searchForLargeStructures>omit large structures</searchForLargeStructures>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>2</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
    <queryId>A263C144</queryId>
    <resultCount>137331</resultCount>
    <runtimeStart>2019-01-16T17:42:58Z</runtimeStart>
    <runtimeMilliseconds>1492</runtimeMilliseconds>
    <containsProtein>Y</containsProtein>
    <containsDna>N</containsDna>
    <containsRna>N</containsRna>
    <containsHybrid>N</containsHybrid>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>3</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.NumberOfChainsQuery</queryType>
    <description>Number of Chains Search : Min Number of Chains=0 Max Number of Chains=1</description>
    <queryId>432D9715</queryId>
    <resultCount>57858</resultCount>
    <runtimeStart>2019-01-16T17:43:00Z</runtimeStart>
    <runtimeMilliseconds>34785</runtimeMilliseconds>
    <struct_asym.numChains.min>0</struct_asym.numChains.min>
    <struct_asym.numChains.max>1</struct_asym.numChains.max>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>4</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.BiolUnitQuery</queryType>
    <description>Oligomeric state Search : Min Number of oligomeric state=0 Max Number of oligomeric state=1</description>
    <queryId>111AF9C9</queryId>
    <resultCount>60014</resultCount>
    <runtimeStart>2019-01-16T17:43:35Z</runtimeStart>
    <runtimeMilliseconds>893</runtimeMilliseconds>
    <oligomeric_statemin>0</oligomeric_statemin>
    <oligomeric_statemax>1</oligomeric_statemax>
  </orgPdbQuery>
 </queryRefinement>
  <queryRefinement>
  <queryRefinementLevel>5</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.TreeQuery</queryType>
    <description>ScopTree Search for Alpha and beta proteins (a/b)</description>
    <queryId>E3E70695</queryId>
    <resultCount>11843</resultCount>
    <runtimeStart>2019-01-16T17:42:10Z</runtimeStart>
    <runtimeMilliseconds>46</runtimeMilliseconds>
    <t>11</t>
    <n>51349</n>
    <nodeDesc>Alpha and beta proteins (a/b)</nodeDesc>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
"""


"""

Query	Structures Found	Conjunction	Structures Displayed	Search time (seconds)
Ligand Search : Has free ligands=no	36116	
36116	0.427
Search for Large Structures : search option = omit large structures	147263	and	35933	1.593
Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid	137331	and	32469	1.492
Number of Chains Search : Min Number of Chains=0 Max Number of Chains=1	57858	and	15149	34.785
Oligomeric state Search : Min Number of oligomeric state=0 Max Number of oligomeric state=1	60014	and	8977	0.893

ScopTree Search for Alpha and beta proteins (a/b)	11843	
11843	0.046
"""




url = 'http://www.rcsb.org/pdb/rest/search'
# data = urllib.quote(queryText).encode("utf-8")
request = urllib.request.Request(url, data=queryText.encode("utf-8"))
response = urllib.request.urlopen(request)
# print()
result = response.read().decode('utf-8')
print(result.count('\n'))
print(result)
