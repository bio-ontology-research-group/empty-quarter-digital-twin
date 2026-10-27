"""Exercise the current process console and public website with real bundles."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright


def main():
    out=Path(__file__).resolve().parents[2]/"reviews/v03/browser"
    out.mkdir(exist_ok=True)
    checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,executable_path="/usr/bin/chromium")
        page=browser.new_page(viewport={"width":1360,"height":1000})
        errors=[];page.on("pageerror",lambda e:errors.append(str(e)))
        page.goto("http://127.0.0.1:8099/")
        page.locator("#sample").wait_for(state="visible")
        page.wait_for_function("!document.querySelector('#sample').disabled")
        page.locator("#run").click();page.get_by_text("Scenario complete.",exact=False).wait_for()
        with page.expect_download() as download:
            page.locator("#download").click()
        downloaded=download.value.path()
        run=json.loads(Path(downloaded).read_text())
        assert run["summary"]["h2_uptake_mol_m2"] > 0
        assert abs(sum(r["salt_balance_residual_mol_m2"] for r in run["trajectory"])) < 1e-12
        checks.append("live_scenario_and_json_download")
        page.screenshot(path=str(out/"console-desktop.png"),full_page=True)
        page.locator("#sample").select_option("F5Dr1")
        assert "Unavailable" in page.locator("#sample-values").inner_text()
        assert "NaN" not in page.locator("#sample-values").inner_text()
        checks.append("missing_metagenomic_observation")
        page.locator("#salt").fill("1");page.locator("#run").click()
        page.get_by_text("Salt exceeded the dilute chemistry domain.",exact=False).wait_for()
        assert "Unavailable" in page.locator("#metrics").inner_text()
        checks.append("high_salinity_abstention")
        page.route("**/v3/simulate",lambda route:route.fulfill(status=503,content_type="application/json",body='{"detail":"Test service unavailable"}'))
        page.locator("#run").click();page.get_by_text("Test service unavailable",exact=False).wait_for()
        assert page.locator("#download").is_disabled()
        assert page.locator("#chart rect").count()==0
        checks.append("failed_run_clears_old_results")
        page.set_viewport_size({"width":390,"height":844})
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        page.screenshot(path=str(out/"console-mobile.png"),full_page=True)
        checks.append("console_mobile_no_overflow")
        page.set_viewport_size({"width":1360,"height":1000})
        page.goto("http://127.0.0.1:5190/model")
        page.locator("#pulse-scenario").wait_for()
        before=page.locator("svg[aria-label^='Computed']").inner_html()
        page.locator("#pulse-scenario").select_option("intact")
        after=page.locator("svg[aria-label^='Computed']").inner_html()
        assert before!=after
        page.locator("#evidence-sample").select_option("F5Dr1")
        assert page.get_by_text("Unavailable",exact=True).count() >= 2
        page.get_by_text("Read all",exact=False).focus();page.keyboard.press("Enter")
        assert page.get_by_text("O1 (ontological).",exact=True).is_visible()
        checks.append("website_scenarios_missingness_and_keyboard_contract")
        page.screenshot(path=str(out/"website-contract.png"),full_page=True)
        page.locator(".twin-v03-page").evaluate("e => { e.scrollTop = 0; }")
        page.screenshot(path=str(out/"website-desktop.png"),full_page=True)
        page.get_by_text("Read all",exact=False).click()
        page.set_viewport_size({"width":390,"height":844})
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        page.locator(".twin-v03-page").evaluate("e => { e.scrollTop = 1000; }")
        assert page.locator(".twin-v03-page").evaluate("e => e.scrollTop > 0")
        page.locator(".twin-v03-page").evaluate("e => { e.scrollTop = 0; }")
        page.screenshot(path=str(out/"website-mobile.png"),full_page=True)
        checks.append("website_mobile_no_overflow")
        page.route("**/twin-data/eq_biodt_v03.json",lambda route:route.fulfill(status=503,body="Unavailable"))
        page.reload();page.get_by_text("Evidence bundle could not be loaded",exact=False).wait_for()
        assert page.locator("#pulse-scenario").count()==0
        checks.append("website_fetch_failure")
        assert not errors,errors
        checks.append("no_browser_javascript_errors")
        browser.close()
    (out/"checks.json").write_text(json.dumps({"passed":checks,"count":len(checks)},indent=2)+"\n")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":main()
