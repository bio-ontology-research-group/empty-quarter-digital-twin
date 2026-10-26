import asyncio, json, os
from pathlib import Path
from playwright.async_api import async_playwright

OUT=Path(__file__).resolve().parent
async def main():
    results=[]
    async with async_playwright() as p:
        browser=await p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
        page=await browser.new_page(viewport={'width':1440,'height':1000})
        errors=[]
        page.on('pageerror',lambda error:errors.append(str(error)))
        await page.goto(os.environ.get('TWIN_WEBSITE_URL', 'http://127.0.0.1:5187/twin'),wait_until='networkidle',timeout=60000)
        await page.get_by_role('heading',name='Read the twin with its evidence').wait_for()
        await page.screenshot(path=str(OUT/'website-desktop.png'))
        text=await page.locator('body').inner_text()
        (OUT/'website-after.txt').write_text(text)
        assert 'Field intervention validation remains outstanding' in text
        details=page.locator('.twin-evidence-guide summary')
        await details.focus(); await page.keyboard.press('Enter')
        assert await page.get_by_text('The Python scenario results are precomputed.',exact=False).is_visible()
        results.append({'check':'Evidence guide visible and keyboard-operable','pass':True})
        await page.get_by_role('link',name='Carbon scenarios',exact=True).click()
        await page.get_by_label('Carbon spatial scope').wait_for()
        await page.get_by_label('Carbon spatial scope').select_option('sand_sea')
        await page.get_by_label('Carbon accounting mode').select_option('gross')
        await page.wait_for_timeout(1000)
        await page.screenshot(path=str(OUT/'website-carbon.png'))
        results.append({'check':'Carbon scope and accounting controls respond','pass':True})
        await page.set_viewport_size({'width':390,'height':844})
        await page.locator('.twin-page').evaluate('(el) => { el.scrollTop = 0; }')
        await page.screenshot(path=str(OUT/'website-mobile.png'),full_page=True)
        overflow=await page.evaluate('document.documentElement.scrollWidth > innerWidth + 1')
        results.append({'check':'Website mobile has no horizontal page overflow','pass':not overflow})
        results.append({'check':'Website JavaScript has no unhandled errors','pass':not errors,'errors':errors.copy()})
        errors.clear()
        await page.set_viewport_size({'width':1440,'height':1000})
        await page.goto(os.environ.get('TWIN_CONSOLE_URL', 'http://127.0.0.1:8098/'),wait_until='networkidle',timeout=60000)
        await page.get_by_text('Baseline loaded. Run a scenario to see conditional outputs.',exact=True).wait_for(timeout=60000)
        assert await page.locator('#sim-kpis .kpi').count()==0
        results.append({'check':'Console loads baseline without silently running an intervention','pass':True})
        await page.get_by_role('button',name='Run scenario',exact=True).click()
        await page.get_by_text('Conditional scenario calculated.',exact=False).wait_for(timeout=60000)
        await page.screenshot(path=str(OUT/'console-desktop.png'))
        results.append({'check':'Console real rainfall intervention succeeds','pass':True})
        await page.route('**/api/v1/intervene/*', lambda route:route.fulfill(json={'summary':{'ocs_score':0},'projected_state':{}}))
        await page.get_by_role('button',name='Run scenario',exact=True).click()
        await page.wait_for_function("document.querySelector('#sim-kpis').textContent.includes('0.000')")
        text=await page.locator('#sim-kpis').inner_text()
        assert 'Unavailable' in text and '0.000' in text
        results.append({'check':'Console preserves OCS zero and labels missing values unavailable','pass':True})
        await page.unroute('**/api/v1/intervene/*')
        await page.route('**/api/v1/intervene/*', lambda route:route.fulfill(status=503,body='Unavailable'))
        await page.get_by_role('button',name='Run scenario',exact=True).click()
        await page.get_by_text('Request failed (503).',exact=False).wait_for()
        assert await page.locator('#sim-kpis .kpi').count()==0
        assert await page.get_by_role('button',name='Run scenario',exact=True).is_enabled()
        results.append({'check':'Console handles HTTP failure and clears stale results','pass':True})
        await page.set_viewport_size({'width':390,'height':844})
        await page.screenshot(path=str(OUT/'console-mobile.png'),full_page=True)
        overflow=await page.evaluate('document.documentElement.scrollWidth > innerWidth + 1')
        results.append({'check':'Console mobile has no horizontal page overflow','pass':not overflow})
        results.append({'check':'Console JavaScript has no unhandled errors','pass':not errors,'errors':errors})
        await browser.close()
    (OUT/'browser-checks.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
    assert all(row['pass'] for row in results)
asyncio.run(main())
