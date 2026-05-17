WTF DOES MY PANEL MEAN?

TLDR: i can read a VO₂ max test in my sleep. then i looked at my own blood panel and realized i had no idea what i was actually looking at. so i went and built the map — more than 130 biomarkers, every one tied to its real hazard ratio, on one comparable scale, fact-checked against the primary literature. then i open-sourced the whole thing. here's what i found — and why your "optimal ranges" are quietly lying to you.

last piece was about the bridge. stay alive, stay metabolically clean, let the next wave of medicine reach you. everything you're doing now — zone 2, sleep, protein — has one job: deliver you in one piece to the medicine that hasn't shipped yet.

fine. but it leaves a question hanging.

how do you actually know if you're metabolically clean?

you get a panel. and then you sit there staring at forty numbers, each with a little flag — "in range," "out of range" — and you feel nothing. no signal. which of these actually matters? which one moves the needle on whether i'm alive and sharp in 2056? the panel won't tell you. your doctor says "looks normal." the influencer says everything is a five-alarm fire.

wtf does my panel actually mean?

I KNOW FITNESS. I DIDN'T KNOW MY PANEL.

honest position: cardiorespiratory fitness is my lane. VO₂ max, ventilatory thresholds, lactate, fat oxidation — i have a mental model. hand me a CPET and i'll tell you a story about that person.

the rest of the panel? apoB, hs-CRP, the metabolic markers, the nine flavors of cholesterol, the alphabet soup of inflammation — i had the exact reaction you have. a black box with a "normal" sticker on it.

so i did for the rest of the panel what i'd already done for fitness: build the mental model. except i didn't want my opinion. i wanted the evidence. so i went and pulled the hazard ratios.

for more than 130 biomarkers, across five outcomes — all-cause mortality, heart disease, cancer, dementia, frailty — what does the actual literature say this number does to your risk? more than 500 evidence cells. every one a hazard ratio with a citation.

then i did the part nobody does. i verified it. every quantified claim, fact-checked against its primary source — twice. when a number couldn't be anchored to a real paper, i deleted it instead of keeping it. dozens of cells died that way. that process alone taught me more than the building did.

and i put the whole thing on github. open source. free. read it, break it, correct it.

THE OPTIMAL RANGE LIE

the first thing that fell apart was the phrase "optimal range."

your lab's reference range is not optimal. it is the middle 95% of everyone who walked in for a blood draw. it means "not unusual." it does not mean "good." normal cholesterol, in a population where heart disease is the leading killer, is not a target — it's the average of a sick crowd.

and it's binary. in-range, out-of-range. a green checkmark or a red flag. but risk isn't a cliff you fall off — it's a ramp you walk up. a number can be flagged "normal" and still be costing you years.

worse, it's a snapshot. one reading tells you a level. it doesn't tell you a direction. and direction is the entire game.

so i threw the ranges out. every marker got re-expressed the honest way: as a hazard ratio per standard deviation — how much your risk actually moves as the number moves. continuous, not binary. comparable across markers measured in completely different units. that's the spine of the whole project.

WHAT THE FIELD OVERSELLS

once everything is on one scale, the hype gets easy to see.

hs-CRP. everyone's favorite inflammation number. strongly predictive — and, by the genetic evidence, basically not causal. lowering CRP itself does nothing; it's a thermometer, not a thermostat. (IL-6, one step upstream, is the real signal. nobody sells you the IL-6 test.)

biological age clocks. the first-generation ones predict mortality weakly. the slick consumer ones vary wildly. GlycanAge has almost no prospective mortality data behind it. telomere length is so measurement-noisy it's nearly useless at the individual level. you are being sold a decimal point of false precision.

micronutrients — the biggest gap in the field between observational hype and reality. vitamin D, selenium, omega-3, the B-vitamins all look fantastic in population studies and then flop in the randomized trials that actually test them. homocysteine predicts beautifully — and lowering it did nothing for events.

and — uncomfortably, in my own lane — VO₂ max. its association with mortality is the strongest of any modifiable marker on the board. but the causal evidence is thinner than the way it gets sold. train anyway. obviously train. just know that part of that signal is fit people being healthier for other reasons.

WHAT TO ACTUALLY WATCH

the short, defensible list:

apoB. the single most defensible number on your panel. it counts the particles that actually cause atherosclerosis, the genetics back causation, it's cheap, and it's deeply modifiable. if you add one test, add this.

blood pressure. causal, free, and most people treat theirs as a footnote.

the metabolic markers — HbA1c, fasting insulin. the engine warning lights.

albuminuria. a cheap urine test almost nobody orders, that quietly predicts cardiovascular death. the most under-rated marker i found.

and your functional markers — fitness, strength, gait speed. the thing i already believed, that the verification confirmed cold: what your body can do predicts more outcomes, more broadly, than almost anything in a vial.

THE POINT

the panel was never the problem. the missing piece is the interpretation layer — the thing that tells you which of the forty numbers is the steering wheel and which are just dashboard lights.

so i built that layer, verified it, and gave it away. go read your own numbers against it.

it won't drive the car for you. but it'll tell you what you're actually looking at.

keep going.
— jack
