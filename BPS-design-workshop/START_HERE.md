## Open the bundle

1. Unzip the bundle.
2. In Codex, open the **BPS-design-workshop** folder itself as the workspace.
   The hidden `.agents` directory must be directly inside the workspace root.
3. Start a new Codex conversation and ask, “What skills are available in this workspace?” Confirm that breeding-scheme-drafter appears. 

If the skill does not appear, confirm that the whole folder was extracted,
then close and reopen the workspace and start a new conversation. 

## Suggested opening prompt

Use `$breeding-scheme-drafter` in design-only mode. R is not installed, and
I do not want implementation or validation during this session. Treat this
as a genuinely new scheme. First have a patient, stage-by-stage conversation
with me to understand the biological process, terminology, tacit decisions,
traits, timing, and which populations must be stored versus can remain
temporary. Help me simplify unnecessary complications. Do not write R code.
After I approve your interpretation, help me agree on the event verbs and
scheduler, then provide the event implementation plan. Stop with a clear
handoff for later implementation. The BPS 0.2.0 source package is available
in `BreedingProgramSimulator/` if selective reference is useful; do not try
to read the entire package before beginning the conversation.

## What design-only mode produces

The conversation should progress through three review points:

1. an agreed interpretation of the breeding process;
2. an agreed event skeleton and scheduler; and
3. an agreed implementation plan describing the main BPS and AlphaSimR
   operations that would later be used.

It deliberately stops before R code, configuration files, package commands,
or smoke tests. The resulting handoff can later be opened on a computer with R
and continued with the same skill in implementation mode.
