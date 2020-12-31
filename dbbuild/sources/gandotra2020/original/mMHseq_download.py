#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import sys
from time import sleep


def initialize(driver, tscale=1.0):
    driver.get("https://mmhseq.shinyapps.io/mMHseq/")
    # Select the appropriate tab, wait for the content to load
    WebDriverWait(driver, 5 * tscale).until(EC.element_to_be_clickable((By.LINK_TEXT, "MH haplotype frequency summary table"))).click()
    WebDriverWait(driver, 5 * tscale).until(EC.presence_of_element_located((By.ID, "mh_frequency_table")))


def parse_marker_elements(driver, tscale=1.0):
    """Click the "marker" dropdown, wait for content to load, and parse list contents"""
    WebDriverWait(driver, 5 * tscale).until(EC.element_to_be_clickable((By.ID, "mh_region4-selectized"))).click()
    markers = WebDriverWait(driver, 5 * tscale).until(EC.presence_of_all_elements_located((By.XPATH, "//div[contains(@class, 'option') and @data-selectable='']")))
    markers = sorted(markers, key=lambda m: m.get_attribute("innerHTML"))
    sleep(3 * tscale)
    markers[-1].click()
    sleep(3 * tscale)
    print(f"[initialize] Found {len(markers)} markers")
    return markers


def parse_population_elements(driver, tscale=1.0):
    """Click the "population" dropdown, wait for content to load, and parse list contents"""
    WebDriverWait(driver, 5 * tscale).until(EC.element_to_be_clickable((By.ID, "population4-selectized"))).click()
    populations = WebDriverWait(driver, 5 * tscale).until(EC.presence_of_all_elements_located((By.XPATH, "//div[contains(@class, 'option') and @data-selectable='']")))
    populations = [p for p in populations if not p.get_attribute("innerHTML").startswith("mh")]
    populations = sorted(populations, key=lambda p: p.get_attribute("innerHTML"))
    sleep(3 * tscale)
    populations[-1].click()
    sleep(3 * tscale)
    print(f"[initialize] Found {len(populations)} populations")
    return populations


def get_dropdowns(driver):
    dropdowns = driver.find_elements_by_xpath("//div[contains(@class, 'selectize-input') and contains(@class, 'has-options') and contains(@class, 'full')]")
    marker_drop = dropdowns[5]
    pop_drop = dropdowns[6]
    return marker_drop, pop_drop


def try_click(elem, tscale=1.0):
    for _ in range(3):
        try:
            elem.click()
            break
        except Exception:
            sleep(1 * tscale)
    else:
        raise RuntimeError(f"Could not click marker {elem.get_attribute('innerHTML')} after 3 attempts")


def main(tscale=1.0, downloaddir=None):
    profile = webdriver.FirefoxProfile()
    if downloaddir:
        ddir = os.path.abspath(downloaddir)
        os.makedirs(ddir, exist_ok=True)
        profile.set_preference("browser.download.dir", ddir)
    profile.set_preference("browser.download.folderList", 2);
    profile.set_preference("browser.download.manager.showWhenStarting", False);
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/csv");
    with webdriver.Firefox(profile) as driver:
        initialize(driver)
        markers = parse_marker_elements(driver)
        populations = parse_population_elements(driver)
        marker_drop, pop_drop = get_dropdowns(driver)
        for marker in markers:
            print(f"[download] Marker {marker.get_attribute('innerHTML')}...", end="", flush=True, file=sys.stderr)
            marker_drop.click()
            try_click(marker)
            for pop in populations:
                pop_drop.click()
                try_click(pop)
                timeout = 5 * tscale
                try:
                    WebDriverWait(driver, timeout).until(EC.element_to_be_clickable((By.LINK_TEXT, "CSV"))).click()
                except Exception:
                    WebDriverWait(driver, timeout).until(EC.element_to_be_clickable((By.LINK_TEXT, "CSV"))).click()
            print("done!", file=sys.stderr)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dir", metavar="D", help="set download directory; by default, Firefox default "
        "directory is used"
    )
    parser.add_argument(
        "-t", "--timeout", metavar="T", type=float, default=1.0, help="multiply the default "
        "network timeout limits by a multiple of T; use '--timeout 2.0' or higher for very slow "
        "connections"
    )
    return parser

if __name__ == "__main__":
    args = cli().parse_args()
    main(tscale=args.timeout, downloaddir=args.dir)
