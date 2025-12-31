import React, { useState, useEffect } from 'react';
import './KeyboardShortcutsHelp.css';
import { KeyboardShortcut } from '../hooks/useKeyboardShortcuts';

interface KeyboardShortcutsHelpProps {
  shortcuts: KeyboardShortcut[];
  isOpen: boolean;
  onClose: () => void;
}

const KeyboardShortcutsHelp: React.FC<KeyboardShortcutsHelpProps> = ({
  shortcuts,
  isOpen,
  onClose,
}) => {
  useEffect(() => {
    if (isOpen) {
      const handleEscape = (e: KeyboardEvent) => {
        if (e.key === 'Escape') {
          onClose();
        }
      };
      window.addEventListener('keydown', handleEscape);
      return () => window.removeEventListener('keydown', handleEscape);
    }
  }, [isOpen, onClose]);

  if (!isOpen) return null;

  const formatKey = (shortcut: KeyboardShortcut): string => {
    const parts: string[] = [];
    if (shortcut.ctrl || shortcut.meta) parts.push('Ctrl');
    if (shortcut.shift) parts.push('Shift');
    if (shortcut.alt) parts.push('Alt');
    parts.push(shortcut.key.toUpperCase());
    return parts.join(' + ');
  };

  return (
    <div className="keyboard-shortcuts-overlay" onClick={onClose}>
      <div className="keyboard-shortcuts-modal" onClick={(e) => e.stopPropagation()}>
        <div className="keyboard-shortcuts-header">
          <h2>⌨️ Keyboard Shortcuts</h2>
          <button className="keyboard-shortcuts-close" onClick={onClose}>
            ✕
          </button>
        </div>
        <div className="keyboard-shortcuts-content">
          {shortcuts.length === 0 ? (
            <p className="keyboard-shortcuts-empty">No shortcuts available</p>
          ) : (
            <div className="keyboard-shortcuts-list">
              {shortcuts.map((shortcut, index) => (
                <div key={index} className="keyboard-shortcut-item">
                  <div className="keyboard-shortcut-keys">
                    {formatKey(shortcut)}
                  </div>
                  <div className="keyboard-shortcut-description">
                    {shortcut.description}
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default KeyboardShortcutsHelp;




