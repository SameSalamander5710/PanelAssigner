# AGENTS.md

## Code and Commit Standards

- Avoid emdash —, unicode arrow → or any unicode characters: ×, … ; use ASCII equivalents instead: -, ->, x, ...
- Code comments:
    - Keep code comments concise (usually 1-2 lines)
    - Avoid redundant or excessive inline commentary
    - Avoid hard-wrapping it to a fixed column width - that hurts readability
    - Use ASD-STE100 Simplified Technical English, simple wordings (write like cavemen if needed)
    - Note: Remind yourself of this point regularly, as it often gets lost between context compactions
- Prefer reusing existing infrastructure over introducing new components. Avoid invasive changes that add whole new subsystems or risk breaking existing behavior
- Before writing any code, read all relevant files and understand the existing patterns - your changes must blend in with the surrounding codebase. 

## Git Commit Rules
- NEVER commit files in the `Data/` folder
- The `Data/` directory is ignored in `.gitignore` for a reason
- Never push commits without explicit approval
