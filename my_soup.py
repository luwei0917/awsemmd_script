from bs4 import BeautifulSoup

with open("domains_summary.cgi") as fp:
    soup = BeautifulSoup(fp, "lxml")
table = soup.find("table", attrs={"tablesorter table_results"})

headings = [th.get_text().strip() for th in table.find("tr").find_all("th")]

print(*headings)

datasets = []
for row in table.find_all("tr")[1:]:
    dataset = (td.get_text().strip() for td in row.find_all("td"))
    print(*dataset)
