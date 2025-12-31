import { useEffect, useCallback } from 'react';

export interface KeyboardShortcut {
  key: string;
  ctrl?: boolean;
  shift?: boolean;
  alt?: boolean;
  meta?: boolean;
  action: () => void;
  description: string;
}

export const useKeyboardShortcuts = (shortcuts: KeyboardShortcut[]) => {
  const handleKeyDown = useCallback(
    (event: KeyboardEvent) => {
      // Ignore if user is typing in an input, textarea, or contenteditable
      const target = event.target as HTMLElement;
      if (
        target.tagName === 'INPUT' ||
        target.tagName === 'TEXTAREA' ||
        target.isContentEditable
      ) {
        return;
      }

      shortcuts.forEach((shortcut) => {
        const ctrlMatch = shortcut.ctrl ? event.ctrlKey || event.metaKey : !event.ctrlKey && !event.metaKey;
        const shiftMatch = shortcut.shift ? event.shiftKey : !event.shiftKey;
        const altMatch = shortcut.alt ? event.altKey : !event.altKey;
        const keyMatch = event.key.toLowerCase() === shortcut.key.toLowerCase();

        if (ctrlMatch && shiftMatch && altMatch && keyMatch) {
          event.preventDefault();
          shortcut.action();
        }
      });
    },
    [shortcuts]
  );

  useEffect(() => {
    window.addEventListener('keydown', handleKeyDown);
    return () => {
      window.removeEventListener('keydown', handleKeyDown);
    };
  }, [handleKeyDown]);
};

// Common keyboard shortcuts for the app
export const createAppShortcuts = (handlers: {
  onSearch?: () => void;
  onExport?: () => void;
  onNextRegion?: () => void;
  onPrevRegion?: () => void;
  onHelp?: () => void;
}): KeyboardShortcut[] => {
  const shortcuts: KeyboardShortcut[] = [];

  if (handlers.onSearch) {
    shortcuts.push({
      key: 'f',
      ctrl: true,
      action: handlers.onSearch,
      description: 'Search regions',
    });
  }

  if (handlers.onExport) {
    shortcuts.push({
      key: 'e',
      ctrl: true,
      action: handlers.onExport,
      description: 'Export data',
    });
  }

  if (handlers.onNextRegion) {
    shortcuts.push({
      key: 'ArrowRight',
      action: handlers.onNextRegion,
      description: 'Next region',
    });
  }

  if (handlers.onPrevRegion) {
    shortcuts.push({
      key: 'ArrowLeft',
      action: handlers.onPrevRegion,
      description: 'Previous region',
    });
  }

  if (handlers.onHelp) {
    shortcuts.push({
      key: '?',
      shift: true,
      action: handlers.onHelp,
      description: 'Show help',
    });
  }

  return shortcuts;
};




