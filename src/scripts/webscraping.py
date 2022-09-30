from bs4 import BeautifulSoup
import requests, lxml, os, json
from collections import defaultdict
import re
import sys

def main():
    search_dict = scrape_one_google_scholar_page()

def scrape_one_google_scholar_page():
    # https://requests.readthedocs.io/en/latest/user/quickstart/#custom-headers
    headers = {
        'User-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36'
    }

    # https://requests.readthedocs.io/en/latest/user/quickstart/#passing-parameters-in-urls
    params = {
        'q': 'ATOX1 "pancreatic progenitors"',  # search query
        'hl': 'en'       # language of the search
    }

    html = requests.get('https://scholar.google.com/scholar', headers=headers, params=params).text
    soup = BeautifulSoup(html, 'lxml')

    # JSON data will be collected here
    data = defaultdict(dict)

    # Container where all needed data is located
    result_count = 0
    for result in soup.select('.gs_r.gs_or.gs_scl'):
        title = result.select_one('.gs_rt').text
        title_link = result.select_one('.gs_rt a')['href']
        publication_info = result.select_one('.gs_a').text
        snippet = result.select_one('.gs_rs').text
        cited_by = result.select_one('#gs_res_ccl_mid .gs_nph+ a')['href']
        try:
            pdf_link = result.select_one('.gs_or_ggsm a:nth-child(1)')['href']
        except: 
            pdf_link = None
        data[title] = {
            'title': title,
            'title_link': title_link,
            'publication_info': publication_info,
            'snippet': snippet,
            'cited_by': f'https://scholar.google.com{cited_by}',
            "pdf_link": pdf_link
        }
        if result_count == 0:
            first_hit = data[title]
        result_count += 1

    print(json.dumps(data, indent = 2, ensure_ascii = False))


    # Get number of matches
    all_text = soup.get_text()
    if re.search("ArticlesScholar[0-9]* results", all_text):
        matching = re.findall("ArticlesScholar[0-9]* results", all_text)
        total_matches = re.sub("ArticlesScholar", "", matching[0])
    elif re.search("ArticlesScholarAbout[0-9]* results", all_text):
        matching = re.findall("ArticlesScholarAbout[0-9]* results", all_text)
        total_matches = re.sub("ArticlesScholarAbout", "", matching[0])
    elif re.search("ArticlesScholar [0-9]* results", all_text):
        matching = re.findall("ArticlesScholar [0-9]* results", all_text)
        total_matches = re.sub("ArticlesScholar ", "", matching[0])
    elif re.search("ArticlesScholarAbout [0-9]* results", all_text):
        matching = re.findall("ArticlesScholarAbout [0-9]* results", all_text)
        total_matches = re.sub("ArticlesScholarAbout ", "", matching[0])
    else:
        print("Unknown wording")


    first_hit["total_matches"] = total_matches

    return(first_hit)


if __name__ == "__main__":
    main()